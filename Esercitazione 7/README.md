## Esercizio 7
Il ciclo di simulazione completa, preceduto dall'equilibrazione, viene eseguito con `./run.sh <phase> <sim_type>` dove `phase` può essere:
- `gas`
- `liquid`
- `solid`

mentre `sim_type` può essere:
- `NVE`
- `NVT`
- `Statistic`
Ciascuno di questi fa riferimento a un file `../INPUT/<sim_type>/input.file` specifico in cui sono elencati i parametri della simulazione vera e propria.

Si può anche eseguire singolarmente il `main` con `./main <input.dat> <y/n>`dove:
- `input.dat` deve essere un file dello stesso tipo di prima
- `y/n` a seconda che si desidera che vengano stampati i file `config.xyz` contenenti le posizioni delle particelle ad ogni step

Il `main` è prodotto con il comando `make`
