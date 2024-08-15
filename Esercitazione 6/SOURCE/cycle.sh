#!/bin/bash


sim_type="$1"
if [[ "$sim_type" = "2" ]]; then
    tot_dir="../METRO"
    input_file="metro.dat"
elif [[ "$sim_type" = "3" ]]; then
    tot_dir="../GIBBS"
    input_file="gibbs.dat"
else
    echo "Tipo di simulazione non valido. Usa '2' o '3'."
    exit 1
fi
# Definisci le liste dei file di input (data) e di output (tot)
data_files=("total_energy.dat" "magnetization.dat" "specific_heat.dat" "susceptibility.dat" "acceptance.dat")
tot_files=("te_tot.dat" "m_tot.dat" "cv_tot.dat" "x_tot.dat" "acc_tot.dat")

input_dir="../INPUT"
output_dir="../OUTPUT"

# Crea la directory se non esiste
mkdir -p "$tot_dir"

# Numero di passi nel ciclo (puoi modificarlo in base alle tue esigenze)
steps="$4"


## Lettura del contenuto del file input.dat in un array senza mapfile
#lines=()
#while IFS= read -r line; do
#    lines+=("$line")
#done < "$input_dir/$input_file"
#
#
## Estrai la prima riga
#line="${lines[0]}"
## Suddividi la riga in colonne usando la tabulazione come delimitatore
#IFS=$'\t' read -r -a columns <<< "$line"
#
## Modifica la seconda colonna con il valore calcolato
#columns[1]="$sim_type"
## Ricompone la riga modificata
#new_line=$(IFS=$'\t'; echo "${columns[*]}")
#lines[0]="$new_line"
## Scrivi il contenuto modificato nel file input.dat
#printf "%s\n" "${lines[@]}" > "$input_dir/$input_file"


# Ciclo per ogni step
for ((i=0; i<steps; i++)); do

    # Calcola il valore da inserire
    value=$(echo "$i * $2 + $3" | bc)
    temp="TEMP                   $value"
    
    if grep -q "TEMP" "$input_dir/$input_file"; then
        sed -i '' "s|TEMP.*|$temp|" "$input_dir/$input_file"
        #echo "Modificata la linea contenente 'TEMP' in $data_file con: '$value'"
    else
        echo "Nessuna linea contenente 'TEMP' trovata in $data_file"
    fi
    
    ./main
    
    for ((j=0; j<${#data_files[@]}; j++)); do
        data_file="$output_dir/${data_files[$j]}"
        output_file="$tot_dir/${tot_files[$j]}"
        if [ ! -f "$output_file" ]; then
            touch "$output_file"
        fi
        # Leggi l'ultima riga del file data.dat
        last_line=$(tail -n 1 "$data_file")

        # Appendi l'ultima riga al file tot.dat preceduta dalla temperatura
        echo -e "$value\t$last_line" >> "$output_file"

        # Opzionale: stampa lo step e la riga appesa per conferma
        echo -e "Step $i: Temperatura inserita: $value - Ultima riga di '${data_files[$j]}' appesa a '${tot_files[$j]}' è [$value\t$last_line]"
    done
done
