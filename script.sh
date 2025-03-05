#!/bin/bash

count=10  # Define o número de vezes que o loop será executado

exe_file="./blde.exe"
seed_base=2


for ((j=40; j<=60; j+=5)); do
    qnt=0
    for c in 20 30 40; do
        ((qnt+=1))
        for p in 5 10 15 20; do

            inFile="instances/instance_${j}_${qnt}_${p}.txt"
            echo "Arquivo de entrada: $inFile"
            outFile="results/result_${j}_${qnt}_${p}.txt"
            echo "Arquivo de saída: $outFile"

            for ((i=1; i<=count; i++)); do
                seed=$((seed_base + i))
                echo "Execução $i com seed=$seed"

                $exe_file -genL 30 -genF 300 -popL 20 -popF 30 -CR 0.7 -F 0.9 -Var 1 -seed $seed -func 20 -inFile "$inFile" -outFile "$outFile"
            done
        done
    done
done

