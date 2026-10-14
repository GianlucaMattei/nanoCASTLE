# nanoCASTLE — come lanciarla

> nanoCASTLE (**nano**pore **C**opy-number **A**nd **ST**ructural variant **L**andscape **E**xplorer) è stato sviluppato con il nome in codice *Nanomatizer*: cartelle, script e file di output dentro l'immagine usano ancora quel nome (`/opt/nanomatizer`, `Nanomatizer_<SAMPLE>_aggregated.xlsx`, ...).

Guida sintetica: i comandi `docker run` e le opzioni. Per i dettagli
(formati depmap, output, note) vedere `TUTORIAL.md`.

**Prima di iniziare**

- Genoma: **solo GRCh38**.
- BAM e genoma devono usare gli stessi nomi di cromosoma (`chr1` in
  entrambi, non `chr1` in uno e `1` nell'altro).
- File che servono: BAM ordinata e FASTA del genoma (gli indici `.bai`
  e `.fai` **non servono**: la pipeline li verifica e li crea
  automaticamente, `TUTORIAL.md` §1) e — **solo se usi `-T`** per
  l'annotazione DepMap — i file DepMap `CRISPRGeneEffect.csv` +
  `Model.csv` (senza `-T` non servono, sezione 5).

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

## L'unica regola

Ogni riga `-v` è `IL_TUO_FILE:PATH_DENTRO_IL_CONTAINER`.

- **Lato sinistro**: il tuo file. Nome e cartella liberi.
- **Lato destro**: copiarlo **esattamente** come scritto.

Sostituire `MioCampione` col nome del campione in tutti i punti
(`SAMPLE=`). Il nome della BAM **dentro** il container è sempre
`BAM.srt.bam` (fisso, non dipende dal campione). Usare lo stesso nome
`SAMPLE` per ogni run dello stesso campione.

---

## 1. Run completa

Tutti i caller (NanoGLADIATOR, Sniffles, cuteSV, Gasoline) +
annotazione. Senza `-T` non ci sono colonne DepMap (sezione 5).

```bash
docker run -it --rm --user $(id -u):$(id -g) \
  -e SAMPLE=MioCampione \
  -v ~/dati/MioCampione.bam:/opt/nanomatizer/Nanomatizer/Test/BAM/BAM.srt.bam \
  -v ~/dati/genoma.fasta:/opt/nanomatizer/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  -v ~/risultati_nanomatizer:/opt/nanomatizer/Nanomatizer/Results \
  -v ~/nanocastle_db:/opt/nanomatizer/userdb \
  gmatteigmattei/nanocastle:1.0.0 \
  bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh
```

Per aggiungere l'annotazione DepMap del tuo tumore: aggiungere
`-T <tumore>` (es. `-T Lung`, `-T "CNS/Brain+Pediatric"` per i soli
modelli pediatrici del lineage, o `-T MBL` per gli ematologici) (i due
CSV DepMap devono essere nella cartella database, sezione 0).

> `Databsese` non è un refuso: la cartella dentro il container si
> chiama davvero così.

I risultati finiscono in `~/risultati_nanomatizer`. Se la run viene
**interrotta**, rilanciare lo stesso comando: i caller già
completati vengono saltati. Se **continua a fallire** allo stesso
step, leggere `STATUS.txt` e `logs/pipeline.log` nella cartella
risultati.

---

## 2. Solo caller (senza annotazione)

Senza file DepMap. I caller si scelgono con `-t` (senza `-t` = tutti
e quattro).

```bash
docker run -it --rm --user $(id -u):$(id -g) \
  -e SAMPLE=MioCampione \
  -v ~/dati/MioCampione.bam:/opt/nanomatizer/Nanomatizer/Test/BAM/BAM.srt.bam \
  -v ~/dati/genoma.fasta:/opt/nanomatizer/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  -v ~/risultati_nanomatizer:/opt/nanomatizer/Nanomatizer/Results \
  -v ~/nanocastle_db:/opt/nanomatizer/userdb \
  gmatteigmattei/nanocastle:1.0.0 \
  bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh -m callers -t sniffles,cutesv
```

---

## 3. Solo annotazione

Lavora sui risultati già in cartella. Senza BAM, senza genoma. Usare
**la stessa cartella risultati** della run caller.

```bash
docker run -it --rm --user $(id -u):$(id -g) \
  -e SAMPLE=MioCampione \
  -v ~/risultati_nanomatizer:/opt/nanomatizer/Nanomatizer/Results \
  -v ~/nanocastle_db:/opt/nanomatizer/userdb \
  gmatteigmattei/nanocastle:1.0.0 \
  bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh -m annotate -t sniffles,cutesv
```

Con `-T <tumore>` (CSV DepMap nella cartella database, sezione 0) l'aggregato
include le colonne DepMap del tumore.

`-t` dice quali risultati tool annotare:

- ogni tool **listato** deve avere i suoi risultati in cartella,
  altrimenti la pipeline si ferma;
- i tool **non listati** restano vuoti nell'xlsx finale, anche se i
  loro risultati sono in cartella.

Senza `-t`, devono essere presenti tutti e 5 gli output: `sniffles`,
`cutesv`, `gladiator`, `gasoline`, `intertrans`.

---

## 4. Aggiungere i caller uno alla volta

Sempre la stessa cartella risultati:

```bash
... -m callers -t sniffles           # comando della sezione 2
... -m callers -t cutesv             # comando della sezione 2
... -m annotate -t sniffles,cutesv   # comando della sezione 3
```

I caller già completati vengono saltati, quindi nulla viene
ricalcolato.

---

## 5. L'annotazione DepMap (opzionale)

**Per default non c'è**: senza `-T` (e senza file custom `-d`)
l'annotazione DepMap è saltata e l'aggregato non ha colonne DepMap.
I due CSV DepMap non servono.

**`-T <tumore>`**: annotazione DepMap sul tumore del campione.
Passare il valore esatto di `OncotreeLineage` del `Model.csv`
(es. `Lung`, `CNS/Brain`, `Ovary/Fallopian Tube`), oppure
`<lineage>+Pediatric` per i soli modelli pediatrici di quel lineage
(es. `CNS/Brain+Pediatric`), oppure `MBL` per gli ematologici.
I due CSV DepMap devono essere nella cartella database (sezione 0):

```bash
  ... run_pipeline_sequential.sh -m annotate -T Lung
```

**File propri** (`-d`, stesso formato di `CRISPRGeneEffect.csv`):
montarli dove si vuole e passare i path dentro il container a `-d`.
I due CSV DepMap standard non servono. `-T` e `-d FILE` non si
possono usare insieme.

```bash
  -v ~/dati/mio_screen.csv:/opt/nanomatizer/depmap/mio_screen.csv \
  ...
  ... run_pipeline_sequential.sh -d /opt/nanomatizer/depmap/mio_screen.csv
```

---

## Opzioni

| Opzione | Valori | Default |
|---|---|---|
| `-m` | `full` (caller + annotazione), `callers`, `annotate` | `full` |
| `-t` | `sniffles`, `cutesv`, `gladiator`, `gasoline`, `intertrans`, `all` (separati da virgola) | `all` |
| `-T` | tipo di tumore (`OncotreeLineage` del Model.csv, `<lineage>+Pediatric`, o `MBL`) — §5 | (nessuno: DepMap saltato) |
| `-d` | `none`, `FILE1,FILE2,...` (file custom) — §5 | `none` |
| `-p` | numero di core da usare | tutti |

- `intertrans` è un secondo output della run Gasoline:
  `-m callers -t intertrans` lancia Gasoline.
- `-T`/`-d` sono ignorati con `-m callers`.

## Output

Nella cartella risultati:

- una cartella per tool;
- CSV/xlsx di annotazione;
- `Nanomatizer_<SAMPLE>_aggregated.xlsx` (tabella finale);
- `STATUS.txt`;
- `logs/`.
