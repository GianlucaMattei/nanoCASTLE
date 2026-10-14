# nanoCASTLE — guida all'uso dell'immagine

> nanoCASTLE (**nano**pore **C**opy-number **A**nd **ST**ructural variant **L**andscape **E**xplorer) è stato sviluppato con il nome in codice *Nanomatizer*: cartelle, script e file di output dentro l'immagine usano ancora quel nome (`/opt/nanomatizer`, `Nanomatizer_<SAMPLE>_aggregated.xlsx`, ...).

Questa guida è per chi **scarica l'immagine Docker** e vuole correre le
analisi sul proprio campione. Non serve installare nulla: tool, pipeline,
database di annotazione e ambienti R sono già nell'immagine — tranne
tre database esterni da scaricare una volta (§0).

**La pipeline è GRCh38-only.** Il genoma deve essere GRCh38 (hg38):
target, GFF, blacklist ENCODE, EnsDb e i nomi dei cromosomi nei caller
sono tutti hg38. Su hg37/GRCh37 la pipeline non funziona.

**Prefisso dei cromosomi**: i nomi dei cromosomi devono essere
**coerenti tra BAM e genoma** — il formato usato in produzione è quello
UCSC/Ensembl con prefisso `chr` (`chr1`, `chr2`, … `chrX`, `chrY`)
sia nei header `@SQ` della BAM sia nell'indice `.fai` del genoma. Se
la tua BAM o il tuo genoma usano nomi senza `chr` (stilo NCBI), la
pipeline non è in grado di trattarli: rigenera i file con i nomi
`chr…`.

Cosa devi avere sul tuo computer (il nome dei file sul tuo disco è
libero, quello che conta è il path **dentro** il container, §1):

- la tua **BAM WGS** (cromosomi primari, GRCh38) — per le run con caller
  (casi 1-4);
- il **genoma GRCh38** `.fa` (es. la no-alt
  `GCA_000001405.15_GRCh38_no_alt_analysis_set` di UCSC/Ensembl) — per
  le run con caller (casi 1-4). L'indice `.fai` **non serve**: la
  pipeline lo verifica e lo crea automaticamente (§1);
- una **cartella risultati** da collegare al container (§1) — sempre
  (non deve necessariamente essere vuota: la pipeline è idempotente,
  vedi caso 1);
- i due CSV **DepMap** (`CRISPRGeneEffect.csv` + `Model.csv`) — solo
  per le run con STEP6 (casi 1, 2, 5, 6) **se passi `-T`** (vedi §4);
  senza `-T` l'annotazione DepMap viene saltata e non servono.

---

## 0. Database esterni (download una tantum)

Alcune risorse **non sono dentro l'immagine** perché i loro termini d'uso non
ne permettono la ridistribuzione (o perché sono grandi e opzionali). Si
scaricano una volta dai provider originali, in una cartella a scelta
(`~/nanocastle_db` negli esempi), e la cartella va montata su
`/opt/nanomatizer/userdb` in ogni run che annota:

| File | Risorsa | Serve per |
|---|---|---|
| `OncoKB.tsv` | OncoKB Cancer Gene List | annotazione (**obbligatorio**) |
| `population_cnv_grch38.txt` | DECIPHER population CNV | annotazione (**obbligatorio**) |
| `trrust_rawdata.human.tsv` | TRRUST v2 | annotazione (**obbligatorio**) |
| `DGV.GS.hg38.gff3` | DGV Gold Standard (hg38, ~380 MB) | annotazione (**obbligatorio**) |
| `CRISPRGeneEffect.csv`, `Model.csv` | DepMap Public | solo con `-T` (opzionale) |

```bash
mkdir -p ~/nanocastle_db
docker run --rm --user $(id -u):$(id -g) \
  -v ~/nanocastle_db:/opt/nanomatizer/userdb \
  gmatteigmattei/nanocastle:1.0.0 bash Nanomatizer/download_databases.sh
```

Aggiungere `--no-depmap` per saltare i due CSV DepMap (~430 MB) se non si
usa `-T`. Per DepMap lo script cerca la release "DepMap Public" più
recente sul portale; se il portale non risponde a uno script (verifica
anti-bot) scarica la 24Q4, l'ultima con URL pubblico stabile. La release
scaricata è scritta in `DEPMAP_RELEASE.txt`. Per usare un'altra release
basta copiare i propri `CRISPRGeneEffect.csv` + `Model.csv` (da
https://depmap.org/portal/data_page/?tab=allData) nella cartella.
Scaricando i file si accettano i termini d'uso di ciascun
provider (vedi `/opt/nanomatizer/THIRD_PARTY.md`). Le run solo caller
(`-m callers`) non hanno bisogno di questa cartella.

---

## 1. I montaggi

I path **dentro** il container sono fissi (gli script li usano così);
il path **sull'host** è libero, basta collegarli con `-v host:container`.

Per le run con caller (casi 1-4) servono:

- la tua BAM al path interno **fisso**
  `/opt/nanomatizer/Nanomatizer/Test/BAM/BAM.srt.bam` — il nome è
  sempre questo, indipendentemente dal campione (§2);
- il tuo genoma GRCh38 al path interno
  `/opt/nanomatizer/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta`;
- la cartella risultati al path interno
  `/opt/nanomatizer/Nanomatizer/Results`.

**Gli indici non sono un'opzione.** La pipeline verifica che
`BAM.srt.bam.bai` e `<genoma>.fai` esistano; se mancano, li crea
automaticamente (con `samtools index` / `samtools faidx`) e riprende.
Se la tua BAM ha già il `.bai`, puoi montarlo accanto (stesso nome +
`.bai`) per saltare l'indice al volo; se non ce l'hai, la pipeline lo
genera da sola. Stesso per il `.fai` del genoma.

Per le run con STEP6 (casi 1, 2, 5, 6) si aggiunge la cartella database
del §0 al path interno `/opt/nanomatizer/userdb` (contiene `OncoKB.tsv`,
`population_cnv_grch38.txt`, `trrust_rawdata.human.tsv`,
`DGV.GS.hg38.gff3` e, per `-T`, `CRISPRGeneEffect.csv` +
`Model.csv`).

I file depmap custom si montano **dove vuoi** dentro il
container e a `-d` passi il path interno (es. monta in
`/opt/nanomatizer/depmap/` e passa
`-d /opt/nanomatizer/depmap/mio1.csv,...`).

Note:

- se i tuoi file si chiamano diversamente, monta comunque col nome
  interno indicato: es.
  `-v ~/dati/prova.bam:/opt/nanomatizer/Nanomatizer/Test/BAM/BAM.srt.bam`.
- i risultati (output, log, `STATUS.txt`) finiscono nella **cartella
  risultati** che monti tu: quando il container finisce, i file sono
  già sulla tua macchina.
- la cartella risultati **non deve essere vuota**: se contiene risultati
  validi di una run precedente, gli step completati vengono saltati
  (idempotenza, caso 1).
- nei casi 5-6 la cartella risultati deve già contenere i risultati dei
  caller prodotti da una run precedente (casi 1-4).

---

## 2. Il nome del campione (`SAMPLE`) — obbligatorio

`-e SAMPLE=<nome>` dice alla pipeline come si chiama il campione. Il
nome entra in **tutti** i path dei risultati (cartelle per-tool, CSV,
nome dell'xlsx finale) e deve essere uguale per ogni caso che corri
sullo stesso esperimento.

Esempio: `-e SAMPLE=MySample`. **Non c'è default**: la pipeline è
generica e si ferma se `SAMPLE` non viene passato.

Il nome della BAM, al contrario, **non** dipende dal campione: dentro
il container è sempre `BAM.srt.bam` (§1).

---

## 3. I sei casi d'uso

La pipeline si lancia così (dentro l'immagine):

```
bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh \
  [-m MODALITA'] [-t TOOLS] [-T TUMOR] [-d DEPMAP] [-p NPROC]
```

- `-m` = `full` (default) | `callers` | `annotate`
- `-t` = `all` (default) | elenco da: `sniffles` | `cutesv` | `gasoline` |
  `intertrans` | `gladiator`
- `-T` = tipo di tumore per l'annotazione DepMap (§4) — **opzionale**;
  senza `-T` (e senza file custom `-d`) l'annotazione DepMap è saltata
- `-d` = `none` (default) | `FILE1,FILE2,...` (file custom, §4)
- `-p` = numero di core da usare (default: tutti i core disponibili)

In `full`/`callers`, `-t` seleziona quali **caller** correre
(`gladiator` = NanoGLADIATOR; `intertrans` implica Gasoline, da cui
deriva). In `annotate`, `-t` seleziona quali **tool** annotare e
aggregare.

### Caso 1 — Run completa

I 4 caller (NanoGLADIATOR → Sniffles → cuteSV → Gasoline) +
annotazione + aggregazione finale. È la run di produzione completa.
Senza `-T` l'aggregato esce **senza** colonne DepMap (§4).

```bash
docker run -it --rm --user $(id -u):$(id -g) \
  -e SAMPLE=MySample \
  -v ~/dati/MySample_no_alt.srt.bam:/opt/nanomatizer/Nanomatizer/Test/BAM/BAM.srt.bam \
  -v ~/dati/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta:/opt/nanomatizer/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  -v ~/risultati_nanomatizer:/opt/nanomatizer/Nanomatizer/Results \
  -v ~/nanocastle_db:/opt/nanomatizer/userdb \
  gmatteigmattei/nanocastle:1.0.0 \
  bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh
```

Per aggiungere l'annotazione DepMap del tuo tumore, aggiungi `-T` (i due
CSV DepMap devono essere nella cartella database, §0), ad esempio per un tumore del polmone:

```bash
  ...
  bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh -T Lung
```

- Tempo: **molte ore** (dipende dalla CPU; cuteSV è il più lento). Per
  un campione nuovo il primo step (NanoGLADIATOR) include la fase
  "prepare" di ~35-40 minuti a 10 processori.
- RAM: lascia libera RAM in abbondanza (cuteSV è il tool più
  esigente; la pipeline richiede ≥ 4 GB liberi prima di partirlo).
- Quando finisce, in `~/risultati_nanomatizer` trovi tutto (§5).
- La run è **idempotente**: se la rilanci (o riparti da una cartella
  risultati già parzialmente compilata), gli step già completati e
  validi vengono saltati e si riparte da dove si era fermato.
- Se un step fallisce la pipeline **si ferma** (non ripulisce nulla):
  leggi `STATUS.txt` e `logs/pipeline.log` nella cartella risultati.
- Lanciare il container **col proprio utente** aggiungendo
  `--user $(id -u):$(id -g)` a `docker run`: i risultati sono così di
  proprietà dell'utente e basta che la cartella risultati sia scrivibile
  da lui. Senza `--user` la pipeline gira come uid/gid 1000 e la cartella
  risultati deve essere scrivibile da quell'uid.

### Caso 2 — Run completa solo su tool selezionabili

Come il caso 1, ma corre **solo i caller scelti** con `-t` e lo STEP6
li annota/aggrega. Esempio: solo Sniffles 2 + NanoGLADIATOR:

```bash
docker run -it --rm --user $(id -u):$(id -g) \
  -e SAMPLE=MySample \
  ... (stessi montaggi del caso 1) ... \
  gmatteigmattei/nanocastle:1.0.0 \
  bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh -m full -t sniffles,gladiator
```

Altri esempi: `-m full -t cutesv`, `-m full -t gasoline,intertrans`.
Gli step caller non selezionati vengono saltati; lo STEP6 produce
l'aggregato con solo i tool selezionati (gli altri slot restano vuoti).

### Caso 3 — Run solo caller (tutti)

Corre i 4 caller e si ferma **prima** dello STEP6: niente annotazione,
niente aggregazione. Non servono i file DepMap.

```bash
docker run -it --rm --user $(id -u):$(id -g) \
  -e SAMPLE=MySample \
  -v ~/dati/MySample_no_alt.srt.bam:/opt/nanomatizer/Nanomatizer/Test/BAM/BAM.srt.bam \
  -v ~/dati/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta:/opt/nanomatizer/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  -v ~/risultati_nanomatizer:/opt/nanomatizer/Nanomatizer/Results \
  -v ~/nanocastle_db:/opt/nanomatizer/userdb \
  gmatteigmattei/nanocastle:1.0.0 \
  bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh -m callers
```

### Caso 4 — Run solo caller selezionabili

Come il caso 3, ma solo i caller scelti con `-t`. Esempio: solo
cuteSV + Gasoline:

```bash
docker run -it --rm --user $(id -u):$(id -g) \
  -e SAMPLE=MySample \
  ... (stessi montaggi del caso 3) ... \
  gmatteigmattei/nanocastle:1.0.0 \
  bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh -m callers -t cutesv,gasoline
```

Se il risultato di un caller esiste già ed è valido, quel caller viene
saltato (idempotenza).

### Caso 5 — Run solo annotazione (tutti i tool)

Nessun caller: prende i risultati **già prodotti** (casi 1-4) nella
cartella risultati e fa lo STEP6 completo (annotazione +
aggregazione) su tutti i tool. Non servono BAM né genoma.

```bash
docker run -it --rm --user $(id -u):$(id -g) \
  -e SAMPLE=MySample \
  -v ~/risultati_nanomatizer:/opt/nanomatizer/Nanomatizer/Results \
  -v ~/nanocastle_db:/opt/nanomatizer/userdb \
  gmatteigmattei/nanocastle:1.0.0 \
  bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh -m annotate
```

Senza `-T` l'aggregato esce senza colonne DepMap. Per annotare il tuo
tumore aggiungi `-T <tumore>` (CSV DepMap nella cartella database, §0).

### Caso 6 — Run solo annotazione selezionabili

Come il caso 5, ma annota/aggrega **solo i tool scelti** con `-t`.
Esempio: solo Sniffles + NanoGLADIATOR:

```bash
docker run -it --rm --user $(id -u):$(id -g) \
  -e SAMPLE=MySample \
  -v ~/risultati_nanomatizer:/opt/nanomatizer/Nanomatizer/Results \
  -v ~/nanocastle_db:/opt/nanomatizer/userdb \
  gmatteigmattei/nanocastle:1.0.0 \
  bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh -m annotate -t sniffles,gladiator
```

I CSV per-tool dei tool esclusi, se presenti, vengono azzerati
(0 righe) prima dell'aggregazione; se assenti, il modulo li sostituisce
con un placeholder vuoto. L'aggregato finale contiene solo i tool
selezionati.

---

## 4. L'annotazione DepMap (`-T` e `-d`)

Lo STEP6 può aggiungere le colonne DepMap all'aggregato. **Per
default non lo fa**: senza `-T` (e senza file custom `-d`)
l'annotazione DepMap viene **saltata** e l'aggregato esce senza
colonne DepMap. Per attivarla scegli una delle due vie (sono
**mutuamente esclusive**):

### 4a. `-T TUMOR` — il tumore del tuo campione

Passa il tipo di tumore:

- il valore esatto di `OncotreeLineage` presente nel `Model.csv` della
  release DepMap che usi (es. `Lung`, `CNS/Brain`,
  `Ovary/Fallopian Tube`, `Ovary`); oppure
- `<lineage>+Pediatric` per il solo sottoinsieme **pediatrico** di quel
  lineage (es. `-T "CNS/Brain+Pediatric"` seleziona i modelli con
  `OncotreeLineage == CNS/Brain` e `AgeCategory == Pediatric`); oppure
- `MBL` per i modelli ematologici (`DepmapModelType == MBL`).

La pipeline legge `Model.csv` + `CRISPRGeneEffect.csv` dalla cartella
database montata su `/opt/nanomatizer/userdb` (§0), seleziona i modelli del tumore e
genera **un gruppo** DepMap con le colonne:

```
DepMap.<TUMOR>.Min  DepMap.<TUMOR>.Q1  DepMap.<TUMOR>.Median
DepMap.<TUMOR>.Q3   DepMap.<TUMOR>.Max  <TUMOR>.Sum.Score  <TUMOR>.Sum.Score.match
```

dove `<TUMOR>` è il nome sanificato (es. `Ovary/Fallopian Tube` →
`Ovary_Fallopian_Tube`). Se il tumore non ha modelli in `Model.csv`, la
pipeline **si ferma** con un messaggio che elenca i lineage disponibili.

Esempio: full run con DepMap su un tumore del polmone:

```bash
docker run -it --rm --user $(id -u):$(id -g) \
  -e SAMPLE=MySample \
  -v ~/dati/MySample_no_alt.srt.bam:/opt/nanomatizer/Nanomatizer/Test/BAM/BAM.srt.bam \
  -v ~/dati/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta:/opt/nanomatizer/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  -v ~/risultati_nanomatizer:/opt/nanomatizer/Nanomatizer/Results \
  -v ~/nanocastle_db:/opt/nanomatizer/userdb \
  gmatteigmattei/nanocastle:1.0.0 \
  bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh -m full -T Lung
```

### 4b. `-d FILE1,FILE2,...` — file custom

Annotazione DepMap con i **tuoi** file, in numero arbitrario (anche 1,
anche 5, quanti ne vuoi). Ogni file ha lo stesso formato di
`CRISPRGeneEffect.csv` (colonna 1 = ModelID, colonne successive = geni
con score). Il path che passi è quello **interno al container** (monta
i file dove vuoi, §1). Per ogni file valido si generano le 5 colonne
`DepMap.<nome>.Min/Q1/Median/Q3/Max`; i file non validi vengono
scartati con un avviso (se tutti sono scartati, l'annotazione DepMap è
vuota, la run prosegue).

Esempio: full con 2 file custom:

```bash
docker run -it --rm --user $(id -u):$(id -g) \
  -e SAMPLE=MySample \
  ... (montaggi BAM/GRCh38/risultati del caso 1) ... \
  -v ~/dati/mio_depmap_1.csv:/opt/nanomatizer/depmap/mio_depmap_1.csv \
  -v ~/dati/mio_depmap_2.csv:/opt/nanomatizer/depmap/mio_depmap_2.csv \
  -v ~/risultati_nanomatizer:/opt/nanomatizer/Nanomatizer/Results \
  -v ~/nanocastle_db:/opt/nanomatizer/userdb \
  gmatteigmattei/nanocastle:1.0.0 \
  bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh -m full \
    -d /opt/nanomatizer/depmap/mio_depmap_1.csv,/opt/nanomatizer/depmap/mio_depmap_2.csv
```

`-T` e `-d FILE` insieme **non** si possono usare: la pipeline si ferma
con un errore.

---

## 5. Output (dove li trovi)

Tutto dentro la cartella risultati che hai montato
(`~/risultati_nanomatizer` nell'esempio; in container
`/opt/nanomatizer/Nanomatizer/Results`):

```
Results/
├── NanoGLADIATOR/<SAMPLE>/R6_W1Mb/      # risultati NanoGLADIATOR (finestre 1 Mb)
├── Sniffles2/<SAMPLE>/                  # risultati Sniffles
├── CuteSV/<SAMPLE>/                     # risultati cuteSV
├── Gasoline/Parse/<SAMPLE>/             # risultati Gasoline
├── <CSV/xlsx di annotazione per tool>   # casi 1, 2, 5, 6
├── Nanomatizer_<SAMPLE>_aggregated.xlsx # aggregato finale (casi 1, 2, 5, 6)
├── STATUS.txt                           # stato: ultimo stop/completamento
└── logs/                                # pipeline.log, step6_v2.log
```

- Se la pipeline è **riuscita**: `STATUS.txt` ha la riga di
  completamento e l'xlsx aggregato è presente (nei casi che lo
  producono).
- Se si è **fermata**: `STATUS.txt` dice perché, `logs/pipeline.log`
  ha il dettaglio. La pulizia di cartelle "sporche" spetta a te: lo
  script non rimuove mai nulla.

---

## 6. Messaggi

Tutti i messaggi di avanzamento e di errore sono scritti in inglese
sul terminale e in `Results/logs/`. `STATUS.txt` registra lo stato
dell'ultimo stop/completamento.

---

## 7. Note

- **GRCh38-only** e **prefisso `chr`** coerente tra BAM e genoma:
  vedi in testa.
- **Finestre di copertura 1 Mb**: la colonna `coverage` dell'aggregato
  e il file `NanoGladiator.W1Mb.<SAMPLE>...exons.csv` usano finestre da
  1 megabase. Se il genoma montato è il GRCh38 standard (la pipeline lo
  riconosce automaticamente dall'indice `.fai`), il target `hg38w1000000`
  già incluso nell'immagine viene usato così com'è; se invece monta una
  variante del genoma con lo stesso assetto cromosomico (stessi 24
  cromosomi `chr`, stesse lunghezze), la pipeline costruisce in autonomia
  le finestre 1 Mb corrispondenti alla prima run (una decina di minuti,
  una sola volta) e le riutilizza nelle run successive. Se il genoma non
  corrisponde a questo assetto, la pipeline si ferma con un messaggio
  esplicito.
- **Indici automatici**: `.bai` e `.fai` vengono verificati e, se
  mancanti, creati dalla pipeline (§1): non sono un'opzione.
- **DepMap opzionale**: senza `-T` (e senza file custom `-d`) le
  colonne DepMap non vengono generate e i due CSV DepMap non servono
  (§4).
- L'utente interno di default è `nanomatizer` (uid 1000); l'immagine
  gira anche con qualunque `--user UID:GID` (consigliato: il proprio,
  vedi caso 1).
- Ogni caller ha la propria verifica strutturale (file presenti e non
  vuoti, header corretti): se una verifica fallisce, la pipeline
  STOPPA e avvisa, senza ripulire nulla.
- La run è idempotente: rilanciarla non rigenera gli step già validi.
- `-p N` limita la run a N core (utile su macchine condivise); il
  default è usare tutti i core disponibili.
