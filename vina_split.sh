#!/bin/bash

# vina_split'in bulunduğu dizin (eğer PATH içindeyse bu satırı silin)
# export PATH=$PATH:/path/to/vina_split

# 1'den 10'a kadar olan klasörler için döngü
for dir in {1..10}; do
    # Eğer klasör varsa içine gir
    if [[ -d "$dir" ]]; then
        # Klasör içindeki tüm .pdbqt dosyaları için döngü
        for file in "$dir"/*.pdbqt; do
            # Dosya varsa, vina_split ile işle
            if [[ -f "$file" ]]; then
                ./vina_split.exe --input "$file"
            fi
        done
    fi
done

echo "İşlem tamamlandı."
