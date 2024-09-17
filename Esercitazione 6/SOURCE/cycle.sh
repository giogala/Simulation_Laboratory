#!/bin/bash
# check point

sim_type="$1"
if [[ "$sim_type" = "2" ]]; then
    tot_dir="../METRO/NO-FIELD"
    input_file="metro.dat"
    echo "Metropolis algorithm"
    equil=false
elif [[ "$sim_type" = "3" ]]; then
    tot_dir="../GIBBS/NO-FIELD"
    input_file="gibbs.dat"
    echo "Gibbs algorithm"
    equil=false
elif [[ "$sim_type" = "2f" ]]; then
    tot_dir="../METRO/FIELD"
    input_file="metro_f.dat"
    echo "Metropolis algorithm"
    equil=false
elif [[ "$sim_type" = "3f" ]]; then
    tot_dir="../GIBBS/FIELD"
    input_file="gibbs_f.dat"
    echo "Gibbs algorithm"
    equil=false
else
    echo -e "Tipo di simulazione non valido. Usa '2' (Metropolis) o '3' (Gibbs) oppure '2f' o '3f' per avere un campo esterno H=0.02.\nUso del programma: ./cycle.sh <sim_type> <Temp_intervall> <T_0> <nof_steps>"
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
    
    mv "$output_dir/CONFIG/config.spin" "$input_dir/CONFIG/config.spin"
    
    for ((j=0; j<${#data_files[@]}; j++)); do
        data_file="$output_dir/${data_files[$j]}"
        output_file="$tot_dir/${tot_files[$j]}"
        
        if $equil; then
            # Leggi l'ultima riga del file data
            last_line=$(tail -n 1 "$data_file")
            # Appendi l'ultima riga al file tot.dat preceduta dalla temperatura
            echo -e "$value\t$last_line" >> "$output_file"
    
            # Stampa lo step e la riga appesa per conferma
            echo -e "Step $i: Temperatura inserita: $value - Ultima riga di '${data_files[$j]}' appesa a '${tot_files[$j]}' \n\t [$value\t$last_line]"
        else
            mv "$data_file" "$tot_dir/${value}_${data_files[$j]}"
            echo -e "Equilibrazione a Temperatura $value \n\t Dati in $tot_dir/${value}_${data_files[$j]}"
        fi
    done
    if ! $equil; then
        # Modifica linea della restart nell'input file
        restart="RESTART                1"
        if grep -q "RESTART" "$input_dir/$input_file"; then
            sed -i '' "s|RESTART.*|$restart|" "$input_dir/$input_file"
            echo "$restart"
        else
            echo "Nessuna linea contenente 'RESTART' trovata in $input_file"
        fi
        equil=$true
        i=$(echo "($i-1)" | bc)
    fi
done

# Modifica linea della restart nell'input file
restart="RESTART                0"
if grep -q "RESTART" "$input_dir/$input_file"; then
    sed -i '' "s|RESTART.*|$restart|" "$input_dir/$input_file"
else
    echo "Nessuna linea contenente 'RESTART' trovata in $input_file"
fi
