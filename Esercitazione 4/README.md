## Esercizio 4
Il ciclo di simulazione completa, preceduto dall'equilibrazione, viene eseguito con `./run.sh <sim_type>` dove il `sim_type` può essere:
- `gas`
- `liquid`
- `solid`

Ciascuno di questi fa riferimento a un file `../INPUT/input.file` specifico in cui sono elencati i parametri della simulazione vera e propria.

Si può anche eseguire singolarmente il `main` con `./main <input.dat> <y/n>`dove:
- `input.dat` deve essere un file dello stesso tipo di prima
- `y/n` a seconda che si desidera che vengano stampati i file `config.xyz` contenenti le posizioni delle particelle ad ogni step

Il `main` è prodotto con il comando `make`
