## Guida rapida all'analisi MSD (Confocale)

Questa guida spiega, in modo semplice e pratico, come usare gli script Python per leggere le traiettorie da un file CSV, calcolare la Mean Square Displacement (MSD) e produrre grafici per l'analisi di dati di microscopia confocale.

Gli script funzionano su file CSV esportati (ad es. da TrackMate 7) con le colonne: Track ID, X position, Y position, Time. Le unità tipiche sono micron (X,Y) e secondi (Time).

Prerequisiti minimi:
- Python 3.9+ (Windows va benissimo)
- Librerie: numpy, pandas, matplotlib (solo per i grafici)

---

## Cosa fanno i file

### 1) `data_reader.py`
- Legge il CSV di input e costruisce, per ogni traccia, un oggetto Trajectory con array NumPy: tempo, x, y.
- Riconosce intestazioni in stile TrackMate (maiuscole/minuscole e varianti comuni) e forza i tipi corretti.
- Ordina i punti per tempo, elimina duplicati temporali e stima il passo temporale Δt di ogni traccia come mediana delle differenze di tempo.
- Funzioni chiave:
	- `read_trajectories_from_csv(path)` → dict {track_id → Trajectory}
	- `estimate_global_time_step(trajectories)` → Δt globale (mediana dei Δt per traccia)

Quando usarlo: sempre, perché tutti gli altri script si appoggiano a questo modulo per leggere i dati in maniera robusta.

#### Come vengono riconosciute e raggruppate le particelle (Track ID)
- Il lettore cerca nel CSV una colonna equivalente a "Track ID" (riconosce anche varianti comuni come "track id", "track_id", "TrackId", ecc.).
- I valori della colonna Track ID vengono convertiti a stringa e ripuliti dagli spazi; questo evita problemi se in alcune righe l’ID è numerico e in altre testuale.
- Tutti i punti con lo stesso Track ID vengono raggruppati (group-by) per formare una singola traiettoria.
- All’interno di ogni gruppo:
	- i punti sono ordinati per tempo crescente;
	- eventuali duplicati di tempo vengono rimossi (si conserva il primo);
	- viene stimato un Δt per traccia come mediana delle differenze temporali.
- L’output è un dizionario: {track_id → Trajectory}. Ogni traiettoria contiene time, x, y, dt e n_points.


---

### 2) `msd_analyzer.py`
- Calcola la MSD in 2D in due modi:
	1. Ensemble-averaged MSD (eaMSD): per ciascun lag n confronta la posizione al tempo n con la posizione iniziale (nessuna media sui punti di partenza), poi media tra le tracce disponibili.
	2. Time-averaged MSD (TAMSD) per **una** singola traccia: media su tutte le finestre possibili all’interno della traccia.
- Unità:
	- τ (lag di tempo) in secondi, usando Δt globale (eaMSD) o della traccia (TAMSD).
	- MSD in unità di posizione al quadrato (es. μm², se X e Y sono in micron).
- Parametro “quanto lag usare”: ora è completamente a scelta dell’utente. Se non si specifica `--max-lag-fraction`, viene usato **tutto l’intervallo disponibile** (fino a N_max − 1).
- Esecuzione come script: stampa Δt, n_max, τ_max e i primi valori τ–MSD (utile per un controllo rapido).

Funzioni principali (utilizzabili anche da altri script):
- `calculate_ensemble_msd(trajectories, max_lag_fraction=None, global_dt=None)`
- `calculate_time_averaged_msd_per_track(track, max_lag_fraction=None, dt_override=None)`

---

#### Struttura interna e funzionamento di `msd_analyzer.py`
`msd_analyzer.py` è un modulo che espone una data-class per i risultati e più funzioni per il calcolo della MSD. La logica è suddivisa in unità riutilizzabili, non in una singola classe monolitica.

Componenti principali:
- `MSDResult` (data-class): raccoglie tau, msd, tracks_per_lag, dt_usato, n_max, numero tracce totali e lunghezza massima.
- Funzioni di supporto:
	- `determine_maximum_lag_steps(N_max, fraction)`: converte la frazione scelta in un numero massimo di passi; limita sempre a 1..N_max-1. Se la frazione è assente o non valida, si usa 1.0 (intervallo completo).
	- `build_tau_array(K, dt)`: costruisce i lag temporali tau = [1..K] * dt.
	- `average_across_trajectories(list_of_arrays)`: media (ignorando NaN) i contributi per lag tra le varie tracce e conta quante tracce hanno contribuito a ciascun lag.

Calcolo dell'eaMSD (ensemble, spostamento rispetto all’inizio):
1) Per ogni traccia e per ciascun lag n, si calcola un singolo spostamento rispetto al punto iniziale: (x[n] - x[0])² + (y[n] - y[0])².
2) I vettori per-traccia (dimensione K) vengono impilati e mediati tra le tracce per ogni lag (nan-safe).
3) I passi n vengono convertiti in tempi con tau = n * Δt, dove Δt è il time-step della traccia.

Calcolo del TAMSD (time-averaged, singola traccia):
- Per ciascun lag n, si calcola la media su tutte le finestre temporali all’interno della traccia: media di (x[i+n]-x[i])² + (y[i+n]-y[i])² per i = 0..N-n-1.
- Il Δt usato è il time-step della traccia
- L’output è un MSDResult relativo a una sola traccia (tracks_per_lag = 1).

Note importanti:
- `max_lag_fraction` è opzionale; se omessa, si usa tutto l’intervallo possibile (N_max-1). Se fornita, deve essere 0 < f ≤ 1 e viene sempre clampata in quel range.
- Le medie ignorano i NaN (tracce troppo corte ai lag alti), così ogni lag usa solo i contributi disponibili.
- Un’intercetta positiva della retta MSD è normale in presenza di rumore di localizzazione o blur di acquisizione.

### 3) `eamsd_plot.py`
- Legge il CSV, calcola la **eaMSD** e salva un grafico PNG su assi lineari (MSD vs τ).
- Mostra anche sul terminale le informazioni principali (Δt, n_max, τ_max, numero di tracce).
- Opzioni:
	- `--max-lag-fraction`: frazione dell’ampiezza della traccia più lunga da usare (0 < f ≤ 1). Se omesso, usa tutto (N_max − 1).
	- `--output`: nome del file immagine di uscita (PNG).

---

### 4) `tamsd_plot.py`
- Calcola e disegna la **TAMSD** per **una** traccia.
- Se non si seleziona una traccia con `--track-id`, usa la prima in ordine.
- Opzioni:
	- `--track-id`: ID della traccia da analizzare (come appare nel CSV).
	- `--max-lag-fraction`: frazione dell’ampiezza della traccia selezionata (0 < f ≤ 1). Se omesso, usa tutto (N − 1).
	- `--output`: nome del file immagine di uscita (PNG).
---

## Come si usano (comandi di esempio)

Di seguito alcuni esempi pronti da copiare in PowerShell (Windows). Sostituisci i percorsi con i tuoi file CSV.

### A) Calcolo rapido dell’eaMSD (solo testo nel terminale)
```powershell
python .\msd_analyzer.py .\sferette240nm_spots.csv
```

per usare solo il 25% dell’intervallo di lag (in punti):
```powershell
python .\msd_analyzer.py .\sferette240nm_spots.csv --max-lag-fraction 0.25
```

Suggerimento: per salvare l’output testuale su file
```powershell
python .\msd_analyzer.py .\sferette240nm_spots.csv > msd_risultati.txt
```

### B) Grafico dell’eaMSD (PNG)
```powershell
python .\eamsd_plot.py .\sferette240nm_spots.csv --output ea_plot.png
```

Con frazione di lag personalizzata (es. 0.3):
```powershell
python .\eamsd_plot.py .\sferette240nm_spots.csv --max-lag-fraction 0.3 --output ea_plot.png
```

### C) Grafico della TAMSD per una traccia
```powershell
python .\tamsd_plot.py .\sferette240nm_spots.csv --output tamsd_plot.png
```

Se conosci l’ID della traccia (ad es. "1"):
```powershell
python .\tamsd_plot.py .\sferette240nm_spots.csv --track-id 1 --output tamsd_track1.png
```

Con frazione di lag personalizzata (es. 0.2):
```powershell
python .\tamsd_plot.py .\sferette240nm_spots.csv --track-id 1 --max-lag-fraction 0.2 --output tamsd_track1.png
```

## Interpretazione rapida dei risultati
- τ (tau): tempo di ritardo (lag) in secondi → asse X.
- MSD: media dei quadrati degli spostamenti in μm² (se le coordinate sono in micron) → asse Y.
- Δt globale: stima robusta del passo temporale comune alle tracce (mediana dei Δt per traccia).
- `--max-lag-fraction`: controlla fin dove spingersi con i lag. Se omesso, usa tutto l’intervallo disponibile. Aumentarlo troppo può ridurre il numero di tracce che contribuiscono alle medie ai lag più lunghi.




