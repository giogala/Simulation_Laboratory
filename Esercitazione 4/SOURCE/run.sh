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

# Estrai la temperatura desiderata
t_fin=$(grep "TEMP" "$input_dir/$input_file" | awk '{print $2}')
echo "$t_fin"

# Crea la directory se non esiste
mkdir -p "$tot_dir"
mkdir -p "$tot_dir/CONFIG"
mkdir -p "$tot_dir/EQUIL"

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
        
        temp="TEMP                   $t_new"
        if grep -q "TEMP" "$input_dir/$input_file"; then
            sed -i '' "s|TEMP.*|$temp|" "$input_dir/$input_file"
        else
            echo "Nessuna linea contenente 'TEMP' trovata in $input_file"
        fi
    fi
    # Sposto in una cartella a parte i dati di equilibrazione
    for ((j=0; j<${#data_files[@]}; j++)); do
        ttin=$(echo "scale=2; $t_in" | bc)
        mv "$output_dir/${data_files[$j]}" "$tot_dir/EQUIL/${ttin}_${data_files[$j]}"
    done
done

# Modifica linea della restart nell'input file
restart="RESTART                1"
if grep -q "RESTART" "$input_dir/$input_file"; then
    sed -i '' "s|RESTART.*|$restart|" "$input_dir/$input_file"
else
    echo "Nessuna linea contenente 'RESTART' trovata in $input_file"
fi
echo "Simulazione avviata"

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
temp="TEMP                   $t_fin"
if grep -q "TEMP" "$input_dir/$input_file"; then
    sed -i '' "s|TEMP.*|$temp|" "$input_dir/$input_file"
else
    echo "Nessuna linea contenente 'TEMP' trovata in $input_file"
fi
restart="RESTART                0"
if grep -q "RESTART" "$input_dir/$input_file"; then
    sed -i '' "s|RESTART.*|$restart|" "$input_dir/$input_file"
else
    echo "Nessuna linea contenente 'RESTART' trovata in $input_file"
fi
