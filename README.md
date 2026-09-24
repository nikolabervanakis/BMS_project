### Estimacija stanja napunjenosti (SoC) baterije pomoću SPKF-a

Studentski projekt u MATLAB-u: implementacija **Sigma-Point Kalmanovog filtra (SPKF)** za procjenu stanja napunjenosti Li-ion ćelije (Samsung 25R, 2.5 Ah) pod dinamičkim opterećenjem.

Algoritam je baziran na skripti iz kolegija ECE5550 profesora Gregoryja Pletta (Lekcija 6.6), prilagođen za model baterije s nelinearnim OCV-om i unutarnjim otporom R_0.



## Rezultati

Filtar je testiran na profilu vožnje koji uključuje gradsku vožnju, stajanje na semaforu, ubrzanje i regenerativno kočenje, uz dodan šum voltmetra.

Namjerno je postavljena kriva početna procjena (stvarna baterija kreće s 85%, a filtar pretpostavlja 50%) kako bi se provjerila konvergencija.

<img width="2748" height="1798" alt="results" src="https://github.com/user-attachments/assets/2286b38c-cfa7-4de0-b73c-fd20a6fb324b" />


* Filtar se ispravi i "zaključa" na pravo stanje unutar prvih 15 sekundi.
* Srednja kvadratna greška (RMSE) u stacionarnom stanju je ispod **0.5%**.
* Prava vrijednost je cijelo vrijeme unutar izračunatih granica pouzdanosti.



## Pokretanje
Pokreni skriptu `bms_spkf_soc.m` u MATLAB-u.
