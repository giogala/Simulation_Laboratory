#!/bin/bash


sim_type="$1"
if [[ "$sim_type" = "2" ]]; then
    tot_dir="../METRO"
    input_file="metro.dat"
    echo "Metropolis algorithm"
    equil=$false
elif [[ "$sim_type" = "3" ]]; then
    tot_dir="../GIBBS"
    input_file="gibbs.dat"
    echo "Gibbs algorithm"
    equil=$false
elif [[ "$sim_type" = "2e" ]]; then
    tot_dir="../METRO_EQUILIBRATION"
    input_file="metro.dat"
    echo "Metropolis Equilibration"
    equil=$true
elif [[ "$sim_type" = "3e" ]]; then
    tot_dir="../GIBBS_EQUILIBRATION"
    input_file="gibbs.dat"
    echo "Gibbs Equilibration"
    equil=$true
else
    echo -e "Tipo di simulazione non valido. Usa '2' (Metropolis) o '3' (Gibbs) oppure '2e' o '3e' per l'equilibrazione.\nUso del programma: ./cycle.sh <sim_type> <Temp_intervall> <T_0> <nof_steps>"
    exit 1
fi
# Liste dei file di singola esecuzione (data) e file totali (tot)
data_files=("total_energy.dat" "magnetization.dat" "specific_heat.dat" "susceptibility.dat" "acceptance.dat")
tot_files=("te_tot.dat" "m_tot.dat" "cv_tot.dat" "x_tot.dat" "acc_tot.dat")

input_dir="../INPUT"
output_dir="../OUTPUT"

# Crea la directory se non esiste
mkdir -p "$tot_dir"

# Numero di passi nel ciclo
steps="$4"

for ((j=0; j<${#tot_files[@]}; j++)); do
    output_file="$tot_dir/${tot_files[$j]}"
    if [ -f "$output_file" ]; then
        rm "$output_file"
    fi
    touch "$output_file"
    echo -e "TEMP\t#BLOCK\tACTUAL\tAVERAGE\tERROR" >> "$output_file"
done
# Ciclo per ogni step
for ((i=1; i<=steps; i++)); do

    # Calcola il valore da inserire
    value=$(echo "($i-1) * $2 + $3" | bc)
    temp="TEMP                   $value"
    if [[ $value == .* ]]; then
        value="0$value"
    fi
    # Modifica linea della temperatura nell'input file
    if grep -q "TEMP" "$input_dir/$input_file"; then
        sed -i '' "s|TEMP.*|$temp|" "$input_dir/$input_file"
    else
        echo "Nessuna linea contenente 'TEMP' trovata in $input_file"
    fi
    
    ./main $input_file
    
    for ((j=0; j<${#data_files[@]}; j++)); do
        data_file="$output_dir/${data_files[$j]}"
        output_file="$tot_dir/${tot_files[$j]}"
        
        if $equil; then
            mv "$data_file" "$tot_dir/${value}_${data_files[$j]}"
            echo -e "Equilibrazione a Temperatura $value \n\t Dati in $tot_dir/${value}_${data_files[$j]}"
        else
            # Leggi l'ultima riga del file data
            last_line=$(tail -n 1 "$data_file")
            # Appendi l'ultima riga al file tot.dat preceduta dalla temperatura
            echo -e "$value\t$last_line" >> "$output_file"
    
            # Stampa lo step e la riga appesa per conferma
            echo -e "Step $i: Temperatura inserita: $value - Ultima riga di '${data_files[$j]}' appesa a '${tot_files[$j]}' \n\t [$value\t$last_line]"
        fi
    done
done
