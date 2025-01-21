#!/bin/bash

# Dosya hash'lerini saklamak için bir associative array (hash table) oluştur
declare -A filehashes

# Dizindeki her dosya için döngü
for file in *; do
    # Sadece normal dosyaları ele al
    if [[ -f "$file" ]]; then
        # Dosyanın MD5 hash değerini hesapla
        hash=$(md5sum "$file" | cut -d ' ' -f 1)
        # Hash değeri zaten mevcutsa, dosya ismini ekle
        filehashes[$hash]="${filehashes[$hash]} $file"
    fi
done

# Hash tablosunda her bir girdi için döngü
for hash in "${!filehashes[@]}"; do
    # Birden fazla dosya aynı hash değerine sahipse, bunları listele
    if [ $(wc -w <<< "${filehashes[$hash]}") -gt 1 ]; then
        echo "Aynı içeriğe sahip dosyalar: ${filehashes[$hash]}"
    fi
done
