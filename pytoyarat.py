# Gerekli modüller
import re

# Kaynak dosya adı
source_file = 'a.py'

# Değiştirilecek dosya isimleri listesi
file_names = [
    "Catechixx", "Chlorogxx", "Coumarix", "Fisetinxx", "ISOQUERxx", 
    "Kaempfexx", "Luteolixxx", "Morin_4_chainb", "Naringexx", 
    "Orobol_chainb", "Scutell", "Taxifol", "Trans-c", "Trans-f", 
    "complex", "complex_su", "hesperetinxx", "ligand", "o-couma", 
    "quercet", "receptor", "rutin", "sinapin", "trans-r"
]

# Her dosya ismi için işlem yap
for name in file_names:
    # Dosya ismini küçük harfe çevir
    lower_name = name.lower()

    # Kaynak dosyayı oku
    with open(source_file, 'r') as file:
        content = file.readlines()

    # İlgili satırlarda dosya isimlerini değiştir
    new_content = []
    for line in content:
        line = line.replace("complex_su.prmtop", f"{lower_name}.prmtop")
        line = line.replace("complex_su.crd", f"{lower_name}.crd")
        line = line.replace("output.nc", f"{lower_name}.nc")
        line = line.replace("checkpoint.chk", f"{lower_name}.chk")
        new_content.append(line)

    # Yeni dosyayı kaydet
    new_file_name = f"{lower_name}.py"
    with open(new_file_name, 'w') as new_file:
        new_file.writelines(new_content)
