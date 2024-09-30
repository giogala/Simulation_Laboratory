#!/bin/bash


sim_type="$1"
if [[ "$sim_type" = "gas" ]]; then
    tot_dir="../OUTPUT/Gas"
    input_file="input.gas"
    echo "Gas System"
    equil="false"
    sigma="0.225"
elif [[ "$sim_type" = "liquid" ]]; then
    tot_dir="../OUTPUT/Liquid"
    input_file="input.liquid"
    echo "Liquid System"
    equil="false"
    sigma="0.0001"
elif [[ "$sim_type" = "solid" ]]; then
    tot_dir="../OUTPUT/Solid"
    input_file="input.solid"
    echo "Solid System"
    equil="false"
    sigma="0.0001"
else
    echo -e "Tipo di simulazione non valido\nUso del programma: ./run.sh <phase> <sim_type>"
    echo -e "phase:\n\t'gas'\n\t'liquid'\n\t'solid'"
    echo -e "sim_type:\n\t'Stat': 500000 blocchi da 1 step\n\t'NVE': NVE simulation\n\t'NVT': NVT simulation"
    exit 1
fi

option="$2"
if [[ "$option" = "Stat" ]]; then
    tot_dir="$tot_dir/Statistic"
    input_dir="../INPUT/Statistic"
elif [[ "$option" = "NVE" ]]; then
    tot_dir="$tot_dir/NVE"
    input_dir="../INPUT/NVE"
elif [[ "$option" = "NVT" ]]; then
    tot_dir="$tot_dir/NVT"
    input_dir="../INPUT/NVT"
else
    input_dir="../INPUT"
fi

data_files=("total_energy.dat" "kinetic_energy.dat" "potential_energy.dat" "temperature.dat" "pressure.dat" "acceptance.dat" "gofr.dat" "output.dat")
output_dir="../OUTPUT"

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
# Raccolgo i parametri da input.dat
t_fin=$(grep "TEMP" "$input_dir/$input_file" | awk '{print $2}') # temperatura estratta da $input.dat
step_zero=$(grep "DELTA" "$input_dir/$input_file" | awk '{print $2}') # parametro ∆ del metropolis estratto da $input.dat
blk=$(grep "NBLOCKS" "$input_dir/$input_file" | awk '{print $2}') # numero di blocchi estratto da $input.dat
spt=$(grep "NSTEPS" "$input_dir/$input_file" | awk '{print $2}') # numero di step per blocco estratto da $input.dat
if [[ "$option" = "NVE" ]]; then
    echo -e "T = $t_fin"
else
    echo -e "∆ = $step_zero"
fi
echo -e "$blk blocchi\n$spt step"

# Crea la directory se non esiste
mkdir -p "$tot_dir"
mkdir -p "$tot_dir/CONFIG"
mkdir -p "$tot_dir/EQUIL"

# Modifico il numero di blocchi e di step per blocco per l'equilibrazione
sub "NBLOCKS                20" "NBLOCKS"
sub "NSTEPS                 2000" "NSTEPS"



if [[ "$option" == "NVE" ]]; then
    while [[ "$equil" == "false" ]]; do
        echo "Equilibrazione"
        ./main "$input_dir" "$input_file" "n"
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
else
    while [[ "$equil" == "false" ]]; do
    
        echo "Equilibrazione"
        ./main "$input_dir" "$input_file" "n"
        echo -e "Eseguita\n"
    
        acc_out=$(tail -n 1 "$output_dir/acceptance.dat" | awk -F'\t' '{print $2}') # accettanza finale della simulazione
        step=$(grep "DELTA" "$input_dir/$input_file" | awk '{print $2}') # parametro ∆ del metropolis estratto da $input.dat
        
        # Calcola lo scostamento (delta) tra l'accettnza finale e 0.50
        delta=$(echo "0.50 - $acc_out" | bc)
    
        # Condizione per verificare l'equilibrio (se delta^2 < sigma^2)
        delta_squared=$(echo "$delta * $delta" | bc)
        
        if (( $(echo "$delta_squared < $sigma" | bc -l) )); then
            equil="true"
            echo -e "Raggiunta accettanza $acc_out con ∆ = $step\nEquilibrazione completata"
            mv "$output_dir/CONFIG/config.xyz" "$input_dir/CONFIG/config_eq.xyz"
            step_new=$step
        else
            # Aggiusta delta in properties.dat
            step_new=$(echo "scale=4; $step * $acc_out / 0.5" | bc) # algoritmo da cambiare
            echo -e "Raggiunta accettanza $acc_out con ∆ = $step\nNuova equilibrazione con ∆ = $step_new"
            
            sub "DELTA                  $step_new" "DELTA"
    
        fi
        # Sposto in una cartella a parte i dati di equilibrazione
        for ((j=0; j<${#data_files[@]}; j++)); do
            ttin=$(echo "scale=2; $step_new" | bc)
            mv "$output_dir/${data_files[$j]}" "$tot_dir/EQUIL/${ttin}_${data_files[$j]}"
        done
    done
fi



# Modifica linea della restart nell'input file
sub "RESTART                1" "RESTART"

# Ristabilisco il numero di blocchi e di step per blocco
sub "NBLOCKS                $blk" "NBLOCKS"
sub "NSTEPS                 $spt" "NSTEPS"


# Esegui veramente la simulazione
echo -e "Simulazione avviata:\n$blk blocchi\n$spt step"
./main "$input_dir" "$input_file" "n"

# Sposto i dati nella cartella apposita
for ((j=0; j<${#data_files[@]}; j++)); do
    mv "$output_dir/${data_files[$j]}" "$tot_dir/"
done

# Sposta i .xyz nella cartella CONFIG apposita
#mv "$output_dir/CONFIG/config_"* "$tot_dir/CONFIG/"

# Sposto posizioni equilibrate nella cartella apposita
mv "$input_dir/CONFIG/config_eq.xyz" "$tot_dir/"

# Resetto il restart in input.dat
sub "RESTART                0" "RESTART"
sub "TEMP                   $t_fin" "TEMP"

