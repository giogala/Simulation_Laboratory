## Esercizio 6
Il ciclo di simulazione completa, con il passaggio da una temperatura all'altra, viene eseguito con `./cycle.sh <sim_type> <Temp_intervall> <T_0> <nof_steps>` dove il `sim_type` può essere:
- `2` o `3` per utilizzare rispettivamente gli algoritmi Metropolis o Gibbs
- `2f` o `3f` allo stesso modo ma in presenza di un campo esterno $H=0.02$

Ciascuno di questi fa riferimento a un file `../INPUT/input.file` specifico in cui sono elencati i parametri della simulazione vera e propria.

Si può anche eseguire singolarmente il `main` con `./main <input.dat>`dove`input.dat` deve essere un file dello stesso tipo di prima


Il `main` è prodotto con il comando `make`
