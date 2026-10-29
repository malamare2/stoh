blocks_per_temp = 11000 - 1000  # jer NbSkip = 1000

with open("f_E.dat") as f:
    lines = f.readlines()

temps = []
current_start = 0
block_count = 0
temp_index = 0

for i, line in enumerate(lines):
    
    # preskoči komentare
    if line.startswith("#"):
        continue

    # prazna linija = kraj temperature
    if line.strip() == "":
        if block_count > 0:
            start = current_start
            end = i - 1
            temps.append((temp_index, start, end, block_count))
            temp_index += 1
            block_count = 0
        continue

    # brojimo blokove
    if block_count == 0:
        current_start = i

    block_count += 1

# zadnja temperatura (ako nema trailing \n\n)
if block_count > 0:
    temps.append((temp_index, current_start, len(lines)-1, block_count))

for t in temps:
    print(f"T{t[0]}: blokovi {t[1]} → {t[2]} | ukupno {t[3]}")