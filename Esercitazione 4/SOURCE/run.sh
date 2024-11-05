#!/bin/bash
# check point

sim_type="$1"
if [[ "$sim_type" = "gas" ]]; then
    tot_dir="../OUTPUT/Gas"
    input_file="input.gas"
    echo "Gas System"
    equil="false"
elif [[ "$sim_type" = "liquid" ]]; then
    tot_dir="../OUTPUT/Liquid"
    input_file="input.liquid"
    echo "Liquid System"
    equil="false"
elif [[ "$sim_type" = "solid" ]]; then
    tot_dir="../OUTPUT/Solid"
    input_file="input.solid"
    echo "Solid System"
    equil="false"
else
    echo -e "Tipo di simulazione non valido. Usa 'gas', 'liquid' o 'solid'.\nUso del programma: ./run.sh <sim_type>"
    exit 1
fi

data_files=("total_energy.dat" "kinetic_energy.dat" "potential_energy.dat" "temperature.dat" "pressure.dat" "acceptance.dat" "output.dat")
input_dir="../INPUT"
output_dir="../OUTPUT"

# Estrai parametri desiderati desiderata
t_fin=$(grep "TEMP" "$input_dir/$input_file" | awk '{print $2}') # temperatura estratta da $input.dat
blk=$(grep "NBLOCKS" "$input_dir/$input_file" | awk '{print $2}') # numero di blocchi estratto da $input.dat
spt=$(grep "NSTEPS" "$input_dir/$input_file" | awk '{print $2}') # numero di step per blocco estratto da $input.dat
echo -e "$t_fin\n$blk\n$spt"

# Definisco una funzione per le sostituzioni in input.dat
if [[ "$OSTYPE" == "darwin"* ]]; then
    SED_CMD="sed -i ''"  # Per macOS
else
    SED_CMD="sed -i"     # Per Linux
fi
sub () {
    str="$1"
    par="$2"
    if grep -q "$par" "$input_dir/$input_file"; then
    $SED_CMD "s|$par.*|$str|" "$input_dir/$input_file"
    else
        echo "Nessuna linea contenente '$par' trovata in $input_file"
    fi
# sed -i '' "s|$par.*|$str|" "$input_dir/$input_file"
}

# Crea la directory se non esiste
mkdir -p "$tot_dir"
mkdir -p "$tot_dir/CONFIG"
mkdir -p "$tot_dir/EQUIL"

# Modifico il numero di blocchi e di step per blocco per l'equilibrazione
sub "NBLOCKS                20" "NBLOCKS"
sub "NSTEPS                 2000" "NSTEPS"

while [[ "$equil" == "false" ]]; do
    echo "Equilibrazione"
    ./main "$input_file" "n"
    echo -e "Eseguita\n"
    t_in=$(grep "TEMPERATURE=" "$output_dir/output.dat" | awk -F'=' '{print $2}' | xargs) # temperatura teorica della simulazione
    t_out=$(tail -n 1 "$output_dir/temperature.dat" | awk -F'\t' '{print $3}') # temperatura reale finale della simulazione
    sigma=$(tail -n 1 "$output_dir/temperature.dat" | awk -F'\t' '{print $4}') # incertezza sulla temperatura
    
    # Calcola lo scostamento (delta) tra la temperatura finale e quella teorica
    delta=$(echo "$t_fin - $t_out" | bc)

    # Condizione per verificare l'equilibrio (se delta^2 < 4 sigma^2)
    delta_squared=$(echo "$delta * $delta" | bc)
    sigma_squared=$(echo "$sigma * $sigma * 4" | bc)
    if (( $(echo "$delta_squared < $sigma_squared" | bc -l) )); then
        equil="true"
        echo -e "Raggiunta T = $t_out ± $sigma. Equilibrazione completata"
        mv "$output_dir/CONFIG/velocities.out" "$input_dir/CONFIG/velocities.in"
        mv "$output_dir/CONFIG/config.xyz" "$input_dir/CONFIG/config_eq.xyz"
    else
        # Aggiusta t_in in properties.dat
        t_new=$(echo "$t_in + $delta" | bc)
        echo "Raggiunta T = $t_out ± $sigma. Nuova equilibrazione a T = $t_new"
        sub "TEMP                   $t_new" "TEMP"
    fi
    # Sposto in una cartella a parte i dati di equilibrazione
    for ((j=0; j<${#data_files[@]}; j++)); do
        ttin=$(echo "scale=2; $t_in" | bc)
        mv "$output_dir/${data_files[$j]}" "$tot_dir/EQUIL/${ttin}_${data_files[$j]}"
    done
done

# Modifica linea della restart nell'input file
sub "RESTART                1" "RESTART"
# Ristabilisco il numero di blocchi e di step per blocco
sub "NBLOCKS                $blk" "NBLOCKS"
sub "NSTEPS                 $spt" "NSTEPS"

echo -e "Simulazione avviata:\n$blk blocchi\n$spt step"

# Esegui veramente la simulazione
./main "$input_file" "y"

# Sposto i dati nella cartella apposita
for ((j=0; j<${#data_files[@]}; j++)); do
    mv "$output_dir/${data_files[$j]}" "$tot_dir/"
done

# Sposta i .xyz nella cartella CONFIG apposita
mv "$output_dir/CONFIG/config_"* "$tot_dir/CONFIG/"

# Sposto posizioni e velocità equilibrate nella cartella apposita
mv "$input_dir/CONFIG/velocities.in" "$tot_dir/"
mv "$input_dir/CONFIG/config_eq.xyz" "$tot_dir/"

# Resetto i valori di temperatura e restart in input.dat
sub "TEMP                   $t_fin" "TEMP"
sub "RESTART                0" "RESTART"
