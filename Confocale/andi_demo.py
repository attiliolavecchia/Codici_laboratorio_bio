import marimo

__generated_with = "0.23.6"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    import numpy as np
    try:
        import plotly.graph_objects as go
        plotly_available = True
        plotly_import_error = ""
    except Exception as exc:
        go = None
        plotly_available = False
        plotly_import_error = str(exc)
    try:
        from andi_datasets.analysis import msd_analysis
        from andi_datasets.datasets_theory import datasets_theory
        from andi_datasets.models_phenom import models_phenom
        andi_available = True
        andi_import_error = ""
    except Exception as exc:
        datasets_theory = None
        models_phenom = None
        msd_analysis = None
        andi_available = False
        andi_import_error = str(exc)

    return (
        andi_available,
        andi_import_error,
        datasets_theory,
        go,
        mo,
        models_phenom,
        msd_analysis,
        plotly_available,
        plotly_import_error,
        np,
    )


@app.cell
def _(mo):
    mo.md(
        r"""
        # AnDi – Interactive Diffusion Demo

        Simulate diffusive trajectories, scrub through time with the frame slider,
        and watch how **taMSD**, **eaMSD** and the **Ergodicity Breaking** parameter 
        change with the model parameters.

        | Symbol | Meaning |
        |--------|---------|
        | $\alpha < 1$ | Subdiffusion |
        | $\alpha = 1$ | Normal Brownian |
        | $\alpha > 1$ | Superdiffusion |
        | EB → 0 | Ergodic|
        | EB finite | Ergodicity broken |
        """
    )
    return


@app.cell
def _(mo):
    model = mo.ui.dropdown(
        options={
            "Brownian Motion": "bm",
            "FBM": "fbm",
            "Lévy Walk": "lw",
            "CTRW": "ctrw",
            "Confinement": "confinement",
            "Immobile Traps": "traps",
        },
        value="bm",
        label="Diffusion model",
        full_width=True,
    )
    n_traj = mo.ui.slider(3, 80, value=15, step=1, label="Trajectories N", show_value=True)
    n_frames = mo.ui.slider(50, 1500, value=300, step=50, label="Simulate T frames", show_value=True)
    box_size = mo.ui.slider(64, 384, value=192, step=16, label="Box size L (px)", show_value=True)
    max_lag_frac = mo.ui.slider(0.1, 0.5, value=0.3, step=0.05, label="Max lag (fraction of T)", show_value=True)
    return box_size, max_lag_frac, model, n_frames, n_traj


@app.cell
def _(mo):
    # All model-specific sliders live here so they are always defined
    # (Marimo requires every variable to exist in exactly one cell)
    r_comp     = mo.ui.slider(2, 20,  value=8,    step=1,    label="Compartment radius r (px)")
    n_comp     = mo.ui.slider(5, 200, value=60,   step=5,    label="Number of compartments Nc", show_value=True)
    trans      = mo.ui.slider(0.0, 1.0, value=0.4, step=0.05, label="Transition probability p")
    d_free     = mo.ui.slider(0.05, 2.0, value=0.5, step=0.05, label="D free")
    d_conf     = mo.ui.slider(0.01, 1.5, value=0.1, step=0.01, label="D confined")
    a_free     = mo.ui.slider(0.5, 2.0, value=1.0,  step=0.05, label="α free")
    a_conf     = mo.ui.slider(0.5, 2.0, value=0.9,  step=0.05, label="α confined")
    alpha_fbm  = mo.ui.slider(0.5, 1.5, value=0.85, step=0.05, label="α  (FBM)")
    alpha_ctrw = mo.ui.slider(0.1, 1.0, value=0.7,  step=0.05, label="α  (CTRW)")
    alpha_lw   = mo.ui.slider(1.0, 2.0, value=1.3,  step=0.05, label="α  (Lévy Walk)")
    d_scale    = mo.ui.slider(0.05, 3.0, value=1.0, step=0.05, label="Diffusion scale D")
    # Immobile Traps sliders
    trap_r   = mo.ui.slider(1, 20,  value=5,    step=1,    label="Trap radius r (px)", show_value=True)
    trap_nt  = mo.ui.slider(5, 150, value=25,   step=5,    label="Number of traps Nt", show_value=True)
    trap_pu  = mo.ui.slider(0.01, 0.5,  value=0.15, step=0.01, label="Unbinding prob Pu")
    trap_pb  = mo.ui.slider(0.001, 0.2, value=0.02, step=0.005, label="Binding prob Pb")
    d_trap   = mo.ui.slider(0.05, 2.0, value=0.5, step=0.05, label="D free (traps)")
    return (
        a_conf, a_free, alpha_ctrw, alpha_fbm, alpha_lw,
        d_conf, d_free, d_scale, n_comp, r_comp, trans,
        trap_r, trap_nt, trap_pu, trap_pb, d_trap,
    )


# ── FRAME SLIDER — separate cell so scrubbing never re-triggers simulation ───
@app.cell
def _(mo):
    current_frame = mo.ui.slider(
        1, 1500, value=50, step=1,
        label="▶  Current frame",
        show_value=True, full_width=True,
    )
    return (current_frame,)


# ── DISPLAY ALL CONTROLS ─────────────────────────────────────────────────────
@app.cell
def _(
    a_conf, a_free, alpha_ctrw, alpha_fbm, alpha_lw,
    andi_available, andi_import_error,
    box_size, current_frame, d_conf, d_free, d_scale,
    d_trap, max_lag_frac, model, mo, n_comp, n_frames, n_traj,
    plotly_available, plotly_import_error,
    r_comp, trans, trap_r, trap_nt, trap_pu, trap_pb,
):
    _info = {
        "bm": mo.callout(
            mo.md("**Brownian Motion** — moto diffusivo classico ($\\alpha=1$), ergodico, EB → 0."),
            kind="neutral",
        ),
        "fbm": mo.callout(
            mo.md("**FBM** — Fractional Brownian Motion; ergodico, $\\alpha<1$: subdiffusione, $\\alpha>1$: superdiffusione."),
            kind="info",
        ),
        "lw": mo.callout(
            mo.md("**Lévy Walk** — salti accoppiati spazio-tempo; $\\alpha>1$ superdiffusione."),
            kind="success",
        ),
        "ctrw": mo.callout(
            mo.md("**CTRW** — Continuous Time Random Walk; non-ergodico (EB finito), $\\alpha<1$ sempre."),
            kind="warn",
        ),
        "confinement": mo.callout(
            mo.md("**Confinamento** — particella salta tra compartimenti circolari con FBM/BM all'interno."),
            kind="info",
        ),
        "traps": mo.callout(
            mo.md("**Trappole immobili** — trappole fisse distribuite casualmente; la particella si lega/slega con prob Pb/Pu."),
            kind="warn",
        ),
    }

    if model.value == "bm":
        _specific = mo.vstack([d_scale], gap=0.3)
    elif model.value == "fbm":
        _specific = mo.vstack([alpha_fbm, d_scale], gap=0.3)
    elif model.value == "lw":
        _specific = mo.vstack([alpha_lw, d_scale], gap=0.3)
    elif model.value == "ctrw":
        _specific = mo.vstack([alpha_ctrw, d_scale], gap=0.3)
    elif model.value == "confinement":
        _specific = mo.vstack([r_comp, n_comp, trans, d_free, d_conf, a_free, a_conf], gap=0.3)
    else:  # traps
        _specific = mo.vstack([trap_r, trap_nt, trap_pu, trap_pb, d_trap], gap=0.3)

    _runtime_notes = []
    if not andi_available:
        _runtime_notes.append(
            mo.callout(
                mo.md(
                    "**WASM fallback attivo**: `andi_datasets` non disponibile nel browser. "
                    "La demo usa una simulazione numerica equivalente in NumPy."
                ),
                kind="warn",
            )
        )
        if model.value in {"confinement", "traps"}:
            _runtime_notes.append(
                mo.callout(
                    mo.md(
                        "Per questo modello, in modalità fallback viene usata un'approssimazione "
                        "Browniana confinata nel box."
                    ),
                    kind="neutral",
                )
            )

    if not plotly_available:
        _runtime_notes.append(
            mo.callout(
                mo.md(
                    "**Plotly non disponibile nel runtime web**: i grafici interattivi sono disattivati."
                ),
                kind="warn",
            )
        )
        if plotly_import_error:
            _runtime_notes.append(
                mo.callout(
                    mo.md(f"Dettaglio import Plotly: `{plotly_import_error}`"),
                    kind="neutral",
                )
            )
        if andi_import_error:
            _runtime_notes.append(
                mo.callout(
                    mo.md(f"Dettaglio import: `{andi_import_error}`"),
                    kind="neutral",
                )
            )

    mo.hstack([
        # Left column: all controls
        mo.vstack([
            mo.md("### 🔬 Model"),
            model,
            _info[model.value],
            *_runtime_notes,
            mo.md("##### Common"),
            n_traj, n_frames, box_size, max_lag_frac,
            mo.md("##### Model parameters"),
            _specific,
            mo.md("---"),
            mo.md("### ⏱ Scrub through time"),
            current_frame,
        ], gap=0.3),
    ], widths=[1])
    return


# ── SIMULATION (does NOT depend on current_frame → no re-run when scrubbing) ─
@app.cell
def _(
    a_conf, a_free, alpha_ctrw, alpha_fbm, alpha_lw,
    andi_available,
    box_size, d_conf, d_free, d_scale, d_trap, datasets_theory,
    max_lag_frac, model, models_phenom, msd_analysis,
    n_comp, n_frames, n_traj, np, r_comp, trans,
    trap_r, trap_nt, trap_pu, trap_pb,
):
    N = int(n_traj.value)
    T = int(n_frames.value)
    L = int(box_size.value)
    comp_centers = None
    comp_radius  = None
    is_bounded   = False   # True for models confined in [0, L]

    def _simulate_fallback(_name, _n, _t, _alpha, _scale, _box):
        _steps = np.random.normal(0.0, 1.0, size=(_n, _t, 2))
        _traj = np.cumsum(_steps, axis=1)

        _tt = np.arange(1, _t + 1, dtype=float)
        _pow = np.power(_tt, (_alpha - 1.0) / 2.0)
        _traj = _traj * _pow[None, :, None] * np.sqrt(max(_scale, 1e-12))

        if _name == "ctrw":
            _move_prob = min(1.0, max(0.05, _alpha))
            _moves = (np.random.random(size=(_n, _t, 1)) < _move_prob).astype(float)
            _traj = np.cumsum(_steps * _moves, axis=1) * np.sqrt(max(_scale, 1e-12))
        elif _name == "lw":
            _bursts = (np.random.random(size=(_n, _t, 1)) < 0.03).astype(float)
            _traj = np.cumsum(_steps * (1.0 + 6.0 * _bursts), axis=1) * np.sqrt(max(_scale, 1e-12))

        _bounded = _name in {"confinement", "traps"}
        if _bounded:
            _traj = _traj - np.min(_traj, axis=(1, 2), keepdims=True)
            _maxv = np.max(_traj, axis=(1, 2), keepdims=True)
            _maxv = np.maximum(_maxv, 1e-9)
            _traj = (_traj / _maxv) * float(_box)

        return _traj, _bounded

    if andi_available and model.value == "confinement":
        is_bounded = True
        comp_radius = int(r_comp.value)
        comp_centers = models_phenom._distribute_circular_compartments(
            Nc=int(n_comp.value), r=comp_radius, L=L,
        )
        _traj_tn2, _labels = models_phenom().confinement(
            N=N, T=T, L=L,
            Ds=[float(d_free.value), float(d_conf.value)],
            alphas=[float(a_free.value), float(a_conf.value)],
            comp_center=comp_centers,
            r=comp_radius,
            trans=float(trans.value),
        )
        traj_nt2 = _traj_tn2.transpose(1, 0, 2)
    elif andi_available and model.value == "traps":
        is_bounded = True
        _traj_tn2, _labels = models_phenom().immobile_traps(
            N=N, T=T, L=L,
            r=int(trap_r.value),
            Pu=float(trap_pu.value),
            Pb=float(trap_pb.value),
            Ds=[float(d_trap.value), 0.0],
            alphas=[1.0, 0.0],
            Nt=int(trap_nt.value),
        )
        traj_nt2 = np.array(_traj_tn2).transpose(1, 0, 2)
    elif andi_available:
        # BM: FBM with alpha=1 (model idx 2); FBM, LW, CTRW: theory models
        _model_idx = {"bm": 2, "fbm": 2, "ctrw": 1, "lw": 3}[model.value]
        _alpha = 1.0 if model.value == "bm" else \
                 float(alpha_fbm.value) if model.value == "fbm" else \
                 float(alpha_ctrw.value) if model.value == "ctrw" else \
                 float(alpha_lw.value)
        _ds = datasets_theory().create_dataset(
            T=T, N_models=N, exponents=[_alpha], models=[_model_idx], dimension=2,
        )
        _x = _ds[:, 2 : 2 + T]
        _y = _ds[:, 2 + T : 2 + 2 * T]
        traj_nt2 = np.stack([_x, _y], axis=-1) * np.sqrt(float(d_scale.value))
    else:
        _alpha = 1.0 if model.value == "bm" else \
                 float(alpha_fbm.value) if model.value == "fbm" else \
                 float(alpha_ctrw.value) if model.value == "ctrw" else \
                 float(alpha_lw.value) if model.value == "lw" else 1.0
        _scale = float(d_scale.value if model.value != "traps" else d_trap.value)
        traj_nt2, is_bounded = _simulate_fallback(model.value, N, T, _alpha, _scale, L)

    # MSD analysis
    _max_lag = max(10, int(T * float(max_lag_frac.value)))
    t_lags = np.arange(1, _max_lag + 1)
    if andi_available:
        _ana = msd_analysis()
        tamsd = _ana.tamsd(traj_nt2, t_lags, dim=2)   # shape (n_lags, N)
        _per_alpha = _ana.get_exponent(traj_nt2, t_lags=t_lags.tolist())
    else:
        tamsd = np.empty((len(t_lags), N), dtype=float)
        for _j, _lag in enumerate(t_lags):
            _disp = traj_nt2[:, _lag:, :] - traj_nt2[:, :-_lag, :]
            tamsd[_j, :] = np.nanmean(np.sum(_disp * _disp, axis=2), axis=1)

        _log_t = np.log(t_lags)
        _per_alpha = np.full(N, np.nan, dtype=float)
        for _i in range(N):
            _mask = np.isfinite(tamsd[:, _i]) & (tamsd[:, _i] > 0)
            if np.count_nonzero(_mask) >= 3:
                _per_alpha[_i] = float(np.polyfit(_log_t[_mask], np.log(tamsd[_mask, _i]), 1)[0])

    tamsd_mean = np.nanmean(tamsd, axis=1)

    eamsd = np.array([
        ((traj_nt2[:, k, :] - traj_nt2[:, 0, :]) ** 2).sum(axis=1).mean()
        for k in t_lags
    ])

    # Per-trajectory alpha
    alpha_ta     = float(np.nanmean(_per_alpha))
    alpha_ta_std = float(np.nanstd(_per_alpha))

    # eaMSD alpha via polyfit on log-log
    _mask_ea = np.isfinite(eamsd) & (eamsd > 0)
    alpha_ea = float(np.polyfit(
        np.log(t_lags[_mask_ea]), np.log(eamsd[_mask_ea]), 1
    )[0]) if _mask_ea.sum() >= 3 else float("nan")

    # Ergodicity Breaking  EB(Δ) = <δ̄²²> / <δ̄²>² − 1
    eb = np.nanmean(tamsd ** 2, axis=1) / (np.nanmean(tamsd, axis=1) ** 2) - 1.0

    return (
        L, N, T,
        alpha_ea, alpha_ta, alpha_ta_std,
        comp_centers, comp_radius, is_bounded,
        eamsd, eb,
        t_lags, tamsd, tamsd_mean,
        traj_nt2,
    )


# ── TRAJECTORY PLOT (fast re-run on frame slider change) ─────────────────────
@app.cell
def _(
    L, N, T,
    alpha_ea, alpha_ta, alpha_ta_std,
    comp_centers, comp_radius, is_bounded,
    plotly_available,
    current_frame, go, model, mo, np, traj_nt2,
):
    if not plotly_available or go is None:
        _out = mo.callout(
            mo.md("Grafico traiettorie non disponibile in questa build web (Plotly mancante)."),
            kind="warn",
        )
    else:
        _frame = min(int(current_frame.value), T)
        _slice = traj_nt2[:, :_frame, :]   # N × frame × 2

        _fig = go.Figure()
        _palette = [f"hsl({int(360 * _i / max(N, 1))},70%,50%)" for _i in range(N)]

        for _i in range(N):
            _xs, _ys = _slice[_i, :, 0], _slice[_i, :, 1]
            _fig.add_trace(go.Scatter(
                x=_xs, y=_ys, mode="lines",
                line={"width": 1.6, "color": _palette[_i]},
                showlegend=False,
            ))
            if len(_xs) > 0:
                _fig.add_trace(go.Scatter(   # start dot
                    x=[_xs[0]], y=[_ys[0]], mode="markers",
                    marker={"symbol": "circle", "size": 6, "color": _palette[_i]},
                    showlegend=False,
                ))
                _fig.add_trace(go.Scatter(   # current position dot
                    x=[_xs[-1]], y=[_ys[-1]], mode="markers",
                    marker={"symbol": "circle", "size": 9, "color": _palette[_i],
                            "line": {"color": "black", "width": 1}},
                    showlegend=False,
                ))

        if model.value == "confinement" and comp_centers is not None:
            for _cx, _cy in comp_centers:
                _fig.add_shape(
                    type="circle", xref="x", yref="y",
                    x0=float(_cx - comp_radius), y0=float(_cy - comp_radius),
                    x1=float(_cx + comp_radius), y1=float(_cy + comp_radius),
                    line={"color": "rgba(25,85,160,0.3)", "width": 1},
                )
            _xr, _yr = [0, L], [0, L]
        elif is_bounded:
            _xr, _yr = [0, L], [0, L]
        else:
            _xmin = float(np.min(traj_nt2[:, :, 0]))
            _xmax = float(np.max(traj_nt2[:, :, 0]))
            _ymin = float(np.min(traj_nt2[:, :, 1]))
            _ymax = float(np.max(traj_nt2[:, :, 1]))
            _px = 0.05 * max(1.0, _xmax - _xmin)
            _py = 0.05 * max(1.0, _ymax - _ymin)
            _xr = [_xmin - _px, _xmax + _px]
            _yr = [_ymin - _py, _ymax + _py]

        _fig.update_layout(
            title=f"{model.value.upper()} | frame {_frame}/{T}",
            xaxis_title="x (px)", yaxis_title="y (px)",
            template="plotly_white", height=520,
            margin={"l": 40, "r": 10, "t": 50, "b": 40},
        )
        _fig.update_xaxes(range=_xr)
        _fig.update_yaxes(range=_yr, scaleanchor="x", scaleratio=1)

        _out = mo.vstack([
            mo.callout(
                mo.md(
                    f"**α (taMSD, mean ± std)** = {alpha_ta:.3f} ± {alpha_ta_std:.3f}"
                    f" &emsp;|&emsp; **α (eaMSD)** = {alpha_ea:.3f}"
                    f" &emsp;|&emsp; frame **{_frame}** / {T}"
                ),
                kind="success" if alpha_ta > 0.9 else "warn" if alpha_ta < 0.6 else "info",
            ),
            _fig,
        ])
    _out
    return


# ── MSD + ERGODICITY BREAKING PLOT ───────────────────────────────────────────
@app.cell
def _(
    alpha_ea, alpha_ta, alpha_ta_std,
    eamsd, eb,
    plotly_available,
    go, model, mo, t_lags, tamsd, tamsd_mean,
):
    if not plotly_available or go is None:
        _out = mo.callout(
            mo.md("Grafici MSD/EB non disponibili in questa build web (Plotly mancante)."),
            kind="warn",
        )
    else:
        # ── Figure 1: taMSD & eaMSD ──────────────────────────────────────────────
        _fig_msd = go.Figure()

        for _i in range(min(tamsd.shape[1], 20)):
            _fig_msd.add_trace(go.Scatter(
                x=t_lags, y=tamsd[:, _i], mode="lines",
                line={"color": "rgba(150,150,220,0.18)", "width": 1},
                showlegend=False,
            ))

        _fig_msd.add_trace(go.Scatter(
            x=t_lags, y=tamsd_mean, mode="lines",
            name=f"taMSD   α = {alpha_ta:.2f} ± {alpha_ta_std:.2f}",
            line={"color": "royalblue", "width": 3},
        ))
        _fig_msd.add_trace(go.Scatter(
            x=t_lags, y=eamsd, mode="lines",
            name=f"eaMSD  α = {alpha_ea:.2f}",
            line={"color": "tomato", "width": 3},
        ))

        _fig_msd.update_xaxes(type="log", title_text="Lag Δ (frames)")
        _fig_msd.update_yaxes(type="log", title_text="MSD")
        _fig_msd.update_layout(
            title=f"taMSD & eaMSD — {model.value.upper()}",
            template="plotly_white", height=400,
            legend={"orientation": "h", "yanchor": "bottom", "y": 1.02},
            margin={"l": 55, "r": 20, "t": 60, "b": 45},
        )

        # ── Figure 2: Ergodicity Breaking ─────────────────────────────────────────
        _fig_eb = go.Figure()

        _fig_eb.add_trace(go.Scatter(
            x=t_lags, y=eb, mode="lines",
            name="EB(Δ)",
            line={"color": "seagreen", "width": 2.5},
        ))
        _fig_eb.add_hline(y=0.0, line_color="silver", line_width=1)
        _eb_bm = 2.0 / (3.0 * tamsd.shape[1])
        _fig_eb.add_hline(
            y=_eb_bm, line_color="sandybrown", line_width=1.5, line_dash="dot",
            annotation_text=f"BM limit = {_eb_bm:.3f}",
            annotation_position="bottom right",
        )

        _fig_eb.update_xaxes(type="log", title_text="Lag Δ (frames)")
        _fig_eb.update_yaxes(title_text="EB(Δ)")
        _fig_eb.update_layout(
            title="Ergodicity Breaking  EB(Δ) = ⟨δ̄²²⟩/⟨δ̄²⟩² − 1",
            template="plotly_white", height=340,
            legend={"orientation": "h", "yanchor": "bottom", "y": 1.02},
            margin={"l": 55, "r": 130, "t": 55, "b": 45},
        )

        _out = mo.vstack([_fig_msd, _fig_eb])
    _out
    return


if __name__ == "__main__":
    app.run()
