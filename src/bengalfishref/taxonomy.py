"""Bangladesh & Bay of Bengal Fish and Commercial Crustacean Checklist Harmonizer (Pillar 1)."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import re
import sqlite3
import time
from typing import Any

import requests

from bengalfishref.schema import CORE_MARKERS, compute_bdtpi, upsert_taxa

# Authoritative curated anchor compendium of flagship Bangladeshi freshwater (IUCN Bangladesh 2015;
# Rahman 2005), Sundarbans estuarine, Bay of Bengal marine fishes, and DoF commercial decapods,
# including historical GenBank synonyms that frequently fragment South Asian eDNA/barcode queries.
# Format: (valid_name, order, family, tax_group, common_en, common_bn, habitat_zones, tier, iucn_bd, synonyms)
CURATED_BANGLADESH_ANCHOR_TAXA: tuple[tuple[str, str, str, str, str, str, str, str, str, list[str]], ...] = (
    # --- Clupeiformes (Hilsa, Shads, Anchovies, Sardines) ---
    ("Tenualosa ilisha", "Clupeiformes", "Dorosomatidae", "Actinopterygii", "Hilsa Shad", "Ilish", "GBM_Freshwater,Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", ["Hilsa ilisha", "Clupea ilisha"]),
    ("Tenualosa toli", "Clupeiformes", "Dorosomatidae", "Actinopterygii", "Toli Shad", "Chandana Ilish", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "VU", ["Hilsa toli"]),
    ("Hilsa kelee", "Clupeiformes", "Dorosomatidae", "Actinopterygii", "Kelee Shad", "Gurja Ilish", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Artisanal_Fishery", "LC", ["Macrura kelee"]),
    ("Gudusia chapra", "Clupeiformes", "Dorosomatidae", "Actinopterygii", "Indian River Shad", "Chapila", "GBM_Freshwater", "Core_Commercial", "VU", ["Clupea chapra"]),
    ("Gonialosa manmina", "Clupeiformes", "Dorosomatidae", "Actinopterygii", "Ganges River Gizzard Shad", "Chamila", "GBM_Freshwater,Sundarbans_Estuarine", "Artisanal_Fishery", "LC", []),
    ("Anodontostoma chacunda", "Clupeiformes", "Dorosomatidae", "Actinopterygii", "Chacunda Gizzard Shad", "Koiputi", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Artisanal_Fishery", "LC", []),
    ("Corica soborna", "Clupeiformes", "Clupeidae", "Actinopterygii", "Ganges River Sprat", "Kachki", "GBM_Freshwater,Sundarbans_Estuarine", "Core_Commercial", "LC", []),
    ("Setipinna phasa", "Clupeiformes", "Engraulidae", "Actinopterygii", "Gangetic Hairfin Anchovy", "Phasa", "GBM_Freshwater,Sundarbans_Estuarine", "Core_Commercial", "LC", []),
    ("Setipinna taty", "Clupeiformes", "Engraulidae", "Actinopterygii", "Scaly Hairfin Anchovy", "Teli Phasa", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Artisanal_Fishery", "LC", []),
    ("Coilia dussumieri", "Clupeiformes", "Engraulidae", "Actinopterygii", "Goldspotted Grenadier Anchovy", "Olua", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Coilia ramcarati", "Clupeiformes", "Engraulidae", "Actinopterygii", "Ramcarat Grenadier Anchovy", "Meghi Olua", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Artisanal_Fishery", "LC", []),
    ("Stolephorus commersonnii", "Clupeiformes", "Engraulidae", "Actinopterygii", "Commerson's Anchovy", "Hada", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Artisanal_Fishery", "LC", []),
    ("Sardinella fimbriata", "Clupeiformes", "Dorosomatidae", "Actinopterygii", "Fringescale Sardinella", "Takhia", "Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Sardinella gibbosa", "Clupeiformes", "Dorosomatidae", "Actinopterygii", "Goldstripe Sardinella", "Champa", "Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Ilisha megaloptera", "Clupeiformes", "Pristigasteridae", "Actinopterygii", "Bigeye Ilisha", "Choukka", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Chirocentrus dorab", "Clupeiformes", "Chirocentridae", "Actinopterygii", "Dorab Wolf-herring", "Karati Chela", "Bay_of_Bengal_Marine", "Artisanal_Fishery", "LC", []),

    # --- Cypriniformes (Major & Minor Carps, Barbs, Danionids, Loaches) ---
    ("Labeo rohita", "Cypriniformes", "Cyprinidae", "Actinopterygii", "Rohu", "Rui", "GBM_Freshwater", "Core_Commercial", "LC", ["Cyprinus rohita"]),
    ("Labeo catla", "Cypriniformes", "Cyprinidae", "Actinopterygii", "Catla", "Katla", "GBM_Freshwater", "Core_Commercial", "LC", ["Catla catla", "Gibelion catla"]),
    ("Cirrhinus mrigala", "Cypriniformes", "Cyprinidae", "Actinopterygii", "Mrigal Carp", "Mrigel", "GBM_Freshwater", "Core_Commercial", "NT", ["Cirrhinus cirrhosus"]),
    ("Labeo calbasu", "Cypriniformes", "Cyprinidae", "Actinopterygii", "Orangefin Labeo", "Kalibaus", "GBM_Freshwater", "Core_Commercial", "LC", []),
    ("Labeo bata", "Cypriniformes", "Cyprinidae", "Actinopterygii", "Bata", "Bata", "GBM_Freshwater", "Core_Commercial", "LC", []),
    ("Labeo gonius", "Cypriniformes", "Cyprinidae", "Actinopterygii", "Kuria Labeo", "Ghonia", "GBM_Freshwater", "Core_Commercial", "NT", []),
    ("Labeo nandina", "Cypriniformes", "Cyprinidae", "Actinopterygii", "Nandi Labeo", "Nandina", "GBM_Freshwater", "Conservation_Endemic", "CR", []),
    ("Labeo pangusia", "Cypriniformes", "Cyprinidae", "Actinopterygii", "Pangusia Labeo", "Ghora Muikha", "GBM_Freshwater", "Conservation_Endemic", "EN", []),
    ("Labeo boga", "Cypriniformes", "Cyprinidae", "Actinopterygii", "Boga Labeo", "Bhangangan", "GBM_Freshwater", "Conservation_Endemic", "CR", []),
    ("Cirrhinus reba", "Cypriniformes", "Cyprinidae", "Actinopterygii", "Reba Carp", "Raikhar", "GBM_Freshwater", "Core_Commercial", "NT", []),
    ("Osteobrama cotio", "Cypriniformes", "Cyprinidae", "Actinopterygii", "Cotio", "Dhela", "GBM_Freshwater", "Artisanal_Fishery", "NT", []),
    ("Amblypharyngodon mola", "Cypriniformes", "Danionidae", "Actinopterygii", "Mola Carplet", "Mola", "GBM_Freshwater", "Core_Commercial", "LC", []),
    ("Esomus danrica", "Cypriniformes", "Danionidae", "Actinopterygii", "Flying Barb", "Darkina", "GBM_Freshwater", "Artisanal_Fishery", "LC", ["Esomus danricus"]),
    ("Rasbora daniconius", "Cypriniformes", "Danionidae", "Actinopterygii", "Slender Rasbora", "Darkina", "GBM_Freshwater", "Artisanal_Fishery", "LC", []),
    ("Salmostoma bacaila", "Cypriniformes", "Danionidae", "Actinopterygii", "Large Razorbelly Minnow", "Katari Chela", "GBM_Freshwater", "Artisanal_Fishery", "LC", ["Oxygaster bacaila"]),
    ("Salmostoma phulo", "Cypriniformes", "Danionidae", "Actinopterygii", "Finescale Razorbelly Minnow", "Fulchela", "GBM_Freshwater", "Artisanal_Fishery", "NT", []),
    ("Cabdio morar", "Cypriniformes", "Danionidae", "Actinopterygii", "Morari", "Piali", "GBM_Freshwater", "Artisanal_Fishery", "VU", ["Aspidoparia morar"]),
    ("Barilius bendelisis", "Cypriniformes", "Danionidae", "Actinopterygii", "Hamilton's Barila", "Joia", "GBM_Freshwater", "Conservation_Endemic", "EN", []),
    ("Danio rerio", "Cypriniformes", "Danionidae", "Actinopterygii", "Zebrafish", "Anju", "GBM_Freshwater", "Bycatch_NonCommercial", "NT", ["Brachydanio rerio"]),
    ("Devario devario", "Cypriniformes", "Danionidae", "Actinopterygii", "Sind Danio", "Banshpati", "GBM_Freshwater", "Artisanal_Fishery", "LC", []),
    ("Puntius sophore", "Cypriniformes", "Cyprinidae", "Actinopterygii", "Pool Barb", "Jat Punti", "GBM_Freshwater", "Core_Commercial", "LC", ["Barbus sophore"]),
    ("Puntius chola", "Cypriniformes", "Cyprinidae", "Actinopterygii", "Swamp Barb", "Chola Punti", "GBM_Freshwater", "Artisanal_Fishery", "LC", []),
    ("Pethia ticto", "Cypriniformes", "Cyprinidae", "Actinopterygii", "Ticto Barb", "Tit Punti", "GBM_Freshwater", "Core_Commercial", "VU", ["Puntius ticto"]),
    ("Pethia conchonius", "Cypriniformes", "Cyprinidae", "Actinopterygii", "Rosy Barb", "Kanchan Punti", "GBM_Freshwater", "Artisanal_Fishery", "LC", ["Puntius conchonius"]),
    ("Pethia gelius", "Cypriniformes", "Cyprinidae", "Actinopterygii", "Golden Dwarf Barb", "Gili Punti", "GBM_Freshwater", "Artisanal_Fishery", "NT", ["Puntius gelius"]),
    ("Pethia phutunio", "Cypriniformes", "Cyprinidae", "Actinopterygii", "Spottedsail Barb", "Phutani Punti", "GBM_Freshwater", "Artisanal_Fishery", "LC", ["Puntius phutunio"]),
    ("Systomus sarana", "Cypriniformes", "Cyprinidae", "Actinopterygii", "Olive Barb", "Sarpunti", "GBM_Freshwater", "Core_Commercial", "NT", ["Puntius sarana"]),
    ("Tor putitora", "Cypriniformes", "Cyprinidae", "Actinopterygii", "Golden Mahseer", "Mohashol", "GBM_Freshwater", "Conservation_Endemic", "EN", []),
    ("Tor tor", "Cypriniformes", "Cyprinidae", "Actinopterygii", "Tor Mahseer", "Mahashol", "GBM_Freshwater", "Conservation_Endemic", "CR", []),
    ("Neolissochilus hexagonolepis", "Cypriniformes", "Cyprinidae", "Actinopterygii", "Copper Mahseer", "Katli", "GBM_Freshwater", "Conservation_Endemic", "EN", []),
    ("Hypophthalmichthys molitrix", "Cypriniformes", "Xenocyprididae", "Actinopterygii", "Silver Carp", "Silver Carp", "GBM_Freshwater", "Core_Commercial", "LC", []),
    ("Hypophthalmichthys nobilis", "Cypriniformes", "Xenocyprididae", "Actinopterygii", "Bighead Carp", "Bighead Carp", "GBM_Freshwater", "Core_Commercial", "LC", ["Aristichthys nobilis"]),
    ("Ctenopharyngodon idella", "Cypriniformes", "Xenocyprididae", "Actinopterygii", "Grass Carp", "Grass Carp", "GBM_Freshwater", "Core_Commercial", "LC", []),
    ("Cyprinus carpio", "Cypriniformes", "Cyprinidae", "Actinopterygii", "Common Carp", "Carpio", "GBM_Freshwater", "Core_Commercial", "LC", []),
    ("Barbonymus gonionotus", "Cypriniformes", "Cyprinidae", "Actinopterygii", "Java Barb", "Thai Sarpunti", "GBM_Freshwater", "Core_Commercial", "LC", ["Puntius gonionotus"]),
    ("Botia dario", "Cypriniformes", "Botiidae", "Actinopterygii", "Bengal Loach", "Rani", "GBM_Freshwater", "Conservation_Endemic", "EN", []),
    ("Botia lohachata", "Cypriniformes", "Botiidae", "Actinopterygii", "Reticulate Loach", "Putul", "GBM_Freshwater", "Conservation_Endemic", "EN", []),
    ("Lepidocephalichthys guntea", "Cypriniformes", "Cobitidae", "Actinopterygii", "Guntea Loach", "Gutum", "GBM_Freshwater", "Core_Commercial", "LC", []),
    ("Canthophrys gongota", "Cypriniformes", "Cobitidae", "Actinopterygii", "Gongota Loach", "Pahari Gutum", "GBM_Freshwater", "Artisanal_Fishery", "NT", ["Somileptes gongota"]),
    ("Paracanthocobitis botia", "Cypriniformes", "Nemacheilidae", "Actinopterygii", "Mottled Loach", "Balichata", "GBM_Freshwater", "Artisanal_Fishery", "LC", ["Acanthocobitis botia"]),

    # --- Osteoglossiformes (Featherbacks) ---
    ("Chitala chitala", "Osteoglossiformes", "Notopteridae", "Actinopterygii", "Clown Knifefish", "Chital", "GBM_Freshwater", "Core_Commercial", "EN", ["Notopterus chitala"]),
    ("Notopterus notopterus", "Osteoglossiformes", "Notopteridae", "Actinopterygii", "Bronze Featherback", "Foli", "GBM_Freshwater", "Core_Commercial", "VU", []),

    # --- Siluriformes (Catfishes: Bagrids, Silurids, Schilbeids, Pangasiids, Clariids, Sisorids) ---
    ("Mystus tengara", "Siluriformes", "Bagridae", "Actinopterygii", "Tengara Catfish", "Bajuri Tengra", "GBM_Freshwater", "Core_Commercial", "LC", []),
    ("Mystus vittatus", "Siluriformes", "Bagridae", "Actinopterygii", "Striped Dwarf Catfish", "Tengra", "GBM_Freshwater,Sundarbans_Estuarine", "Core_Commercial", "LC", []),
    ("Mystus cavasius", "Siluriformes", "Bagridae", "Actinopterygii", "Gangetic Mystus", "Golsha Tengra", "GBM_Freshwater,Sundarbans_Estuarine", "Core_Commercial", "NT", []),
    ("Mystus gulio", "Siluriformes", "Bagridae", "Actinopterygii", "Long Whiskers Catfish", "Nona Tengra", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "NT", []),
    ("Mystus bleekeri", "Siluriformes", "Bagridae", "Actinopterygii", "Day's Mystus", "Golsha", "GBM_Freshwater", "Artisanal_Fishery", "LC", []),
    ("Sperata aor", "Siluriformes", "Bagridae", "Actinopterygii", "Long-whiskered Catfish", "Ayre", "GBM_Freshwater", "Core_Commercial", "VU", ["Mystus aor", "Aorichthys aor"]),
    ("Sperata seenghala", "Siluriformes", "Bagridae", "Actinopterygii", "Giant River-catfish", "Guizza Ayre", "GBM_Freshwater", "Core_Commercial", "VU", ["Aorichthys seenghala"]),
    ("Rita rita", "Siluriformes", "Bagridae", "Actinopterygii", "Rita", "Rita", "GBM_Freshwater,Sundarbans_Estuarine", "Core_Commercial", "EN", []),
    ("Wallago attu", "Siluriformes", "Siluridae", "Actinopterygii", "Helicopter Catfish", "Boal", "GBM_Freshwater", "Core_Commercial", "VU", []),
    ("Ompok pabda", "Siluriformes", "Siluridae", "Actinopterygii", "Pabdah Catfish", "Madhu Pabda", "GBM_Freshwater", "Core_Commercial", "EN", []),
    ("Ompok bimaculatus", "Siluriformes", "Siluridae", "Actinopterygii", "Butter Catfish", "Kani Pabda", "GBM_Freshwater", "Core_Commercial", "EN", []),
    ("Ompok pabo", "Siluriformes", "Siluridae", "Actinopterygii", "Pabo Catfish", "Pabo", "GBM_Freshwater", "Conservation_Endemic", "CR", []),
    ("Ailia coila", "Siluriformes", "Ailiidae", "Actinopterygii", "Gangetic Ailia", "Kajuli", "GBM_Freshwater", "Core_Commercial", "LC", []),
    ("Clupisoma garua", "Siluriformes", "Ailiidae", "Actinopterygii", "Garua Bachcha", "Ghaura", "GBM_Freshwater", "Core_Commercial", "EN", []),
    ("Eutropiichthys vacha", "Siluriformes", "Ailiidae", "Actinopterygii", "Batchwa Vacha", "Bacha", "GBM_Freshwater", "Core_Commercial", "LC", []),
    ("Silonia silondia", "Siluriformes", "Ailiidae", "Actinopterygii", "Silond Catfish", "Shilong", "GBM_Freshwater,Sundarbans_Estuarine", "Conservation_Endemic", "EN", []),
    ("Pachypterus atherinoides", "Siluriformes", "Horabagridae", "Actinopterygii", "Indian Potasi", "Batasi", "GBM_Freshwater", "Core_Commercial", "LC", ["Neotropius atherinoides", "Pseudeutropius atherinoides"]),
    ("Pangasius pangasius", "Siluriformes", "Pangasiidae", "Actinopterygii", "Yellowtail Catfish", "Deshi Pangas", "GBM_Freshwater,Sundarbans_Estuarine", "Core_Commercial", "EN", []),
    ("Pangasianodon hypophthalmus", "Siluriformes", "Pangasiidae", "Actinopterygii", "Striped Catfish", "Thai Pangas", "GBM_Freshwater", "Core_Commercial", "LC", ["Pangasius hypophthalmus", "Pangasius sutchi"]),
    ("Clarias magur", "Siluriformes", "Clariidae", "Actinopterygii", "Walking Catfish", "Magur", "GBM_Freshwater", "Core_Commercial", "LC", ["Clarias batrachus"]),
    ("Clarias gariepinus", "Siluriformes", "Clariidae", "Actinopterygii", "North African Catfish", "African Magur", "GBM_Freshwater", "Artisanal_Fishery", "LC", []),
    ("Heteropneustes fossilis", "Siluriformes", "Heteropneustidae", "Actinopterygii", "Stinging Catfish", "Shing", "GBM_Freshwater", "Core_Commercial", "LC", ["Saccobranchus fossilis"]),
    ("Bagarius bagarius", "Siluriformes", "Sisoridae", "Actinopterygii", "Dwarf Goonch", "Baghair", "GBM_Freshwater", "Conservation_Endemic", "CR", []),
    ("Gagata cenia", "Siluriformes", "Sisoridae", "Actinopterygii", "Indian Gagata", "Cenia", "GBM_Freshwater", "Artisanal_Fishery", "LC", []),
    ("Chaca chaca", "Siluriformes", "Chacidae", "Actinopterygii", "Squarehead Catfish", "Cheka", "GBM_Freshwater", "Conservation_Endemic", "EN", []),
    ("Plotosus canius", "Siluriformes", "Plotosidae", "Actinopterygii", "Gray Eel-catfish", "Gang Magur", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "NT", []),
    ("Arius arius", "Siluriformes", "Ariidae", "Actinopterygii", "Threadfin Sea Catfish", "Kata", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Osteogeneiosus militaris", "Siluriformes", "Ariidae", "Actinopterygii", "Soldier Catfish", "Apuia", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Artisanal_Fishery", "LC", []),

    # --- Anabantiformes, Synbranchiformes, Cichliformes & Gobiiformes ---
    ("Anabas testudineus", "Anabantiformes", "Anabantidae", "Actinopterygii", "Climbing Perch", "Koi", "GBM_Freshwater", "Core_Commercial", "LC", []),
    ("Trichogaster fasciata", "Anabantiformes", "Osphronemidae", "Actinopterygii", "Banded Gourami", "Khalisha", "GBM_Freshwater", "Core_Commercial", "LC", ["Colisa fasciata", "Colisa fasciatus"]),
    ("Trichogaster lalius", "Anabantiformes", "Osphronemidae", "Actinopterygii", "Dwarf Gourami", "Lal Khalisha", "GBM_Freshwater", "Artisanal_Fishery", "LC", ["Colisa lalia"]),
    ("Channa striata", "Anabantiformes", "Channidae", "Actinopterygii", "Striped Snakehead", "Shol", "GBM_Freshwater", "Core_Commercial", "LC", ["Ophicephalus striatus"]),
    ("Channa marulius", "Anabantiformes", "Channidae", "Actinopterygii", "Great Snakehead", "Gajar", "GBM_Freshwater", "Core_Commercial", "EN", []),
    ("Channa punctata", "Anabantiformes", "Channidae", "Actinopterygii", "Spotted Snakehead", "Taki", "GBM_Freshwater", "Core_Commercial", "LC", ["Channa punctatus"]),
    ("Channa orientalis", "Anabantiformes", "Channidae", "Actinopterygii", "Walking Snakehead", "Cheng", "GBM_Freshwater", "Artisanal_Fishery", "VU", ["Channa gachua"]),
    ("Channa barca", "Anabantiformes", "Channidae", "Actinopterygii", "Barca Snakehead", "Pipla Shol", "GBM_Freshwater", "Conservation_Endemic", "CR", []),
    ("Nandus nandus", "Anabantiformes", "Nandidae", "Actinopterygii", "Gangetic Leaffish", "Meni / Bheda", "GBM_Freshwater", "Core_Commercial", "NT", []),
    ("Badis badis", "Anabantiformes", "Badidae", "Actinopterygii", "Blue Perch", "Napit Koi", "GBM_Freshwater", "Artisanal_Fishery", "NT", []),
    ("Monopterus cuchia", "Synbranchiformes", "Synbranchidae", "Actinopterygii", "Gangetic Mud Eel", "Kuchia", "GBM_Freshwater,Sundarbans_Estuarine", "Core_Commercial", "VU", ["Ophisternon cuchia"]),
    ("Macrognathus aculeatus", "Synbranchiformes", "Mastacembelidae", "Actinopterygii", "Lesser Spiny Eel", "Tara Baim", "GBM_Freshwater", "Core_Commercial", "NT", []),
    ("Macrognathus pancalus", "Synbranchiformes", "Mastacembelidae", "Actinopterygii", "Barred Spiny Eel", "Guchi Baim", "GBM_Freshwater", "Core_Commercial", "LC", ["Mastacembelus pancalus"]),
    ("Mastacembelus armatus", "Synbranchiformes", "Mastacembelidae", "Actinopterygii", "Zig-zag Eel", "Shal Baim", "GBM_Freshwater", "Core_Commercial", "EN", []),
    ("Oreochromis niloticus", "Cichliformes", "Cichlidae", "Actinopterygii", "Nile Tilapia", "Nilotica / Tilapia", "GBM_Freshwater,Sundarbans_Estuarine", "Core_Commercial", "LC", ["Tilapia nilotica"]),
    ("Oreochromis mossambicus", "Cichliformes", "Cichlidae", "Actinopterygii", "Mozambique Tilapia", "Tilapia", "GBM_Freshwater,Sundarbans_Estuarine", "Core_Commercial", "LC", []),
    ("Glossogobius giuris", "Gobiiformes", "Gobiidae", "Actinopterygii", "Tank Goby", "Bele", "GBM_Freshwater,Sundarbans_Estuarine", "Core_Commercial", "LC", []),
    ("Apocryptes bato", "Gobiiformes", "Gobiidae", "Actinopterygii", "Chewa", "Chiring", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Artisanal_Fishery", "LC", []),
    ("Pseudapocryptes elongatus", "Gobiiformes", "Gobiidae", "Actinopterygii", "Lanceolate Goby", "Chewa", "Sundarbans_Estuarine", "Core_Commercial", "LC", []),
    ("Boleophthalmus boddarti", "Gobiiformes", "Gobiidae", "Actinopterygii", "Boddart's Goggle-eyed Goby", "Dahuk", "Sundarbans_Estuarine", "Artisanal_Fishery", "LC", []),
    ("Odontamblyopus rubicundus", "Gobiiformes", "Gobiidae", "Actinopterygii", "Rubicundus Eelgoby", "Lal Chewa", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Artisanal_Fishery", "LC", []),

    # --- Estuarine & Bay of Bengal Marine Finfishes (Perches, Pomfrets, Croakers, Mackerels, Threadfins) ---
    ("Lates calcarifer", "Carangiformes", "Latidae", "Actinopterygii", "Barramundi / Asian Seabass", "Koral / Bhetki", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Harpadon nehereus", "Aulopiformes", "Synodontidae", "Actinopterygii", "Bombay-duck", "Loitta", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "NT", []),
    ("Pampus argenteus", "Scombriformes", "Stromateidae", "Actinopterygii", "Silver Pomfret", "Foli Chanda / Rupchanda", "Bay_of_Bengal_Marine", "Core_Commercial", "VU", []),
    ("Pampus chinensis", "Scombriformes", "Stromateidae", "Actinopterygii", "Chinese Silver Pomfret", "Rupchanda", "Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Parastromateus niger", "Carangiformes", "Carangidae", "Actinopterygii", "Black Pomfret", "Kalo Chanda", "Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Rastrelliger kanagurta", "Scombriformes", "Scombridae", "Actinopterygii", "Indian Mackerel", "Champa", "Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Scomberomorus commerson", "Scombriformes", "Scombridae", "Actinopterygii", "Narrow-barred Spanish Mackerel", "Maitta", "Bay_of_Bengal_Marine", "Core_Commercial", "NT", []),
    ("Scomberomorus guttatus", "Scombriformes", "Scombridae", "Actinopterygii", "Indo-Pacific King Mackerel", "Bijram", "Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Euthynnus affinis", "Scombriformes", "Scombridae", "Actinopterygii", "Kawakawa", "Bom Maitta", "Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Katsuwonus pelamis", "Scombriformes", "Scombridae", "Actinopterygii", "Skipjack Tuna", "Tuna", "Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Thunnus tonggol", "Scombriformes", "Scombridae", "Actinopterygii", "Longtail Tuna", "Lomba Tuna", "Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Trichiurus lepturus", "Scombriformes", "Trichiuridae", "Actinopterygii", "Largehead Hairtail", "Chhuri", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Lepturacanthus savala", "Scombriformes", "Trichiuridae", "Actinopterygii", "Savalani Hairtail", "Rupali Chhuri", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Polynemus paradiseus", "Carangiformes", "Polynemidae", "Actinopterygii", "Paradise Threadfin", "Tapasi", "GBM_Freshwater,Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Eleutheronema tetradactylum", "Carangiformes", "Polynemidae", "Actinopterygii", "Fourfinger Threadfin", "Tailla", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "VU", []),
    ("Leptomelanosoma indicum", "Carangiformes", "Polynemidae", "Actinopterygii", "Indian Threadfin", "Lakhua", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "VU", ["Polydactylus indicus"]),
    ("Sillaginopsis panijus", "Eupercaria", "Sillaginidae", "Actinopterygii", "Flathead Sillago", "Tular Dandi", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Sillago sihama", "Eupercaria", "Sillaginidae", "Actinopterygii", "Silver Sillago", "Hundra", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Artisanal_Fishery", "LC", []),
    ("Scatophagus argus", "Eupercaria", "Scatophagidae", "Actinopterygii", "Spotted Scat", "Bishtara", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Terapon jarbua", "Centrarchiformes", "Terapontidae", "Actinopterygii", "Jarbua Terapon", "Borgi", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Artisanal_Fishery", "LC", []),
    ("Mugil cephalus", "Mugiliformes", "Mugilidae", "Actinopterygii", "Flathead Grey Mullet", "Bhangon", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Planiliza subviridis", "Mugiliformes", "Mugilidae", "Actinopterygii", "Greenback Mullet", "Bata", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", ["Chelon subviridis", "Liza subviridis"]),
    ("Planiliza parsia", "Mugiliformes", "Mugilidae", "Actinopterygii", "Goldspot Mullet", "Parshe", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", ["Liza parsia"]),
    ("minimugil cascasia", "Mugiliformes", "Mugilidae", "Actinopterygii", "Yellowtail Mullet", "Kachki Bata", "GBM_Freshwater", "Artisanal_Fishery", "VU", ["Sicamugil cascasia"]),
    ("Rhinomugil corsula", "Mugiliformes", "Mugilidae", "Actinopterygii", "Corsula", "Khorsula", "GBM_Freshwater,Sundarbans_Estuarine", "Core_Commercial", "LC", []),
    ("Johnius belangerii", "Eupercaria", "Sciaenidae", "Actinopterygii", "Belanger's Croaker", "Poa", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Otolithoides pama", "Eupercaria", "Sciaenidae", "Actinopterygii", "Pama Croaker", "Lomba Poa", "GBM_Freshwater,Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Protonibea diacanthus", "Eupercaria", "Sciaenidae", "Actinopterygii", "Blackspotted Croaker", "Kala Datina", "Bay_of_Bengal_Marine", "Core_Commercial", "NT", []),
    ("Lutjanus johnii", "Eupercaria", "Lutjanidae", "Actinopterygii", "John's Snapper", "Ranga Choukka", "Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Lutjanus argentimaculatus", "Eupercaria", "Lutjanidae", "Actinopterygii", "Mangrove Red Snapper", "Lal Snapper", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Epinephelus coioides", "Perciformes", "Serranidae", "Actinopterygii", "Orange-spotted Grouper", "Bolia", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Nemipterus japonicus", "Eupercaria", "Nemipteridae", "Actinopterygii", "Japanese Threadfin Bream", "Rupban", "Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Upeneus sulphureus", "Syngnathiformes", "Mullidae", "Actinopterygii", "Sulphur Goatfish", "Sonali Bata", "Bay_of_Bengal_Marine", "Artisanal_Fishery", "LC", []),
    ("Megalops cyprinoides", "Elopiformes", "Megalopidae", "Actinopterygii", "Indo-Pacific Tarpon", "Nanchil", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Artisanal_Fishery", "LC", []),
    ("Chanda nama", "Ovalentaria", "Ambassidae", "Actinopterygii", "Elongate Glass-perchlet", "Lomba Chanda", "GBM_Freshwater", "Core_Commercial", "LC", []),
    ("Parambassis ranga", "Ovalentaria", "Ambassidae", "Actinopterygii", "Indian Glassy Fish", "Ranga Chanda", "GBM_Freshwater", "Core_Commercial", "LC", []),
    ("Xenentodon cancila", "Beloniformes", "Belonidae", "Actinopterygii", "Freshwater Garfish", "Kakila", "GBM_Freshwater", "Core_Commercial", "LC", []),
    ("Hyporhamphus limbatus", "Beloniformes", "Hemiramphidae", "Actinopterygii", "Congaturi Halfbeak", "Ekthota", "GBM_Freshwater,Sundarbans_Estuarine", "Artisanal_Fishery", "LC", []),

    # --- Elasmobranchii (Bay of Bengal & Estuarine Sharks, Guitarfishes, Stingrays) ---
    ("Scoliodon laticaudus", "Carcharhiniformes", "Carcharhinidae", "Elasmobranchii", "Spadenose Shark", "Thuti Hangor", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Artisanal_Fishery", "NT", []),
    ("Carcharhinus limbatus", "Carcharhiniformes", "Carcharhinidae", "Elasmobranchii", "Blacktip Shark", "Kalo Hangor", "Bay_of_Bengal_Marine", "Artisanal_Fishery", "VU", []),
    ("Carcharhinus sorrah", "Carcharhiniformes", "Carcharhinidae", "Elasmobranchii", "Spot-tail Shark", "Hangor", "Bay_of_Bengal_Marine", "Artisanal_Fishery", "NT", []),
    ("Glyphis gangeticus", "Carcharhiniformes", "Carcharhinidae", "Elasmobranchii", "Ganges Shark", "Gang Hangor", "GBM_Freshwater,Sundarbans_Estuarine", "Conservation_Endemic", "CR", []),
    ("Sphyrna lewini", "Carcharhiniformes", "Sphyrnidae", "Elasmobranchii", "Scalloped Hammerhead", "Haturi Hangor", "Bay_of_Bengal_Marine", "Conservation_Endemic", "CR", []),
    ("Glaucostegus granulatus", "Rhinopristiformes", "Glaucostegidae", "Elasmobranchii", "Sharpnose Guitarfish", "Pitambori", "Bay_of_Bengal_Marine", "Conservation_Endemic", "CR", ["Rhinobatos granulatus"]),
    ("Pristis pristis", "Rhinopristiformes", "Pristidae", "Elasmobranchii", "Largetooth Sawfish", "Karat Hangor", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Conservation_Endemic", "CR", ["Pristis microdon"]),
    ("Brevitrygon imbricata", "Myliobatiformes", "Dasyatidae", "Elasmobranchii", "Scaly Whipray", "Shapla Pata", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Artisanal_Fishery", "VU", ["Himantura imbricata"]),
    ("Fluvitrygon fluviatilis", "Myliobatiformes", "Dasyatidae", "Elasmobranchii", "Gangetic Stingray", "Shapla Pata", "GBM_Freshwater,Sundarbans_Estuarine", "Conservation_Endemic", "EN", ["Pastinachus sephen"]),
    ("Maculabatis gerrardi", "Myliobatiformes", "Dasyatidae", "Elasmobranchii", "Sharpnose Stingray", "Shapla Pata", "Bay_of_Bengal_Marine", "Artisanal_Fishery", "EN", ["Himantura gerrardi"]),

    # --- Commercial Crustaceans (Penaeid Shrimps, Palaemonid Prawns, Portunid Crabs, Lobsters) ---
    ("Penaeus monodon", "Decapoda", "Penaeidae", "Crustacea", "Giant Tiger Prawn", "Bagda Chingri", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Penaeus indicus", "Decapoda", "Penaeidae", "Crustacea", "Indian White Prawn", "Chaka Chingri", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", ["Fenneropenaeus indicus"]),
    ("Penaeus merguiensis", "Decapoda", "Penaeidae", "Crustacea", "Banana Prawn", "BagaChingri", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", ["Fenneropenaeus merguiensis"]),
    ("Penaeus semisulcatus", "Decapoda", "Penaeidae", "Crustacea", "Green Tiger Prawn", "Baghtara Chingri", "Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Penaeus japonicus", "Decapoda", "Penaeidae", "Crustacea", "Kuruma Prawn", "Dorakata Chingri", "Bay_of_Bengal_Marine", "Core_Commercial", "LC", ["Marsupenaeus japonicus"]),
    ("Penaeus vannamei", "Decapoda", "Penaeidae", "Crustacea", "Whiteleg Shrimp", "Vannamei Chingri", "Sundarbans_Estuarine", "Core_Commercial", "LC", ["Litopenaeus vannamei"]),
    ("Metapenaeus monoceros", "Decapoda", "Penaeidae", "Crustacea", "Speckled Shrimp", "Horina Chingri", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Metapenaeus brevicornis", "Decapoda", "Penaeidae", "Crustacea", "Yellow Shrimp", "Honne Chingri", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Metapenaeus dobsoni", "Decapoda", "Penaeidae", "Crustacea", "Kadal Shrimp", "Thengua Chingri", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Metapenaeus affinis", "Decapoda", "Penaeidae", "Crustacea", "Jinga Shrimp", "Karaney Chingri", "Bay_of_Bengal_Marine", "Artisanal_Fishery", "LC", []),
    ("Parapenaeopsis sculptilis", "Decapoda", "Penaeidae", "Crustacea", "Rainbow Shrimp", "Ruda Chingri", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Parapenaeopsis stylifera", "Decapoda", "Penaeidae", "Crustacea", "Kiddi Shrimp", "Khudi Chingri", "Bay_of_Bengal_Marine", "Artisanal_Fishery", "LC", []),
    ("Acetes indicus", "Decapoda", "Sergestidae", "Crustacea", "Jawla Paste Shrimp", "Gura Ichha", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Macrobrachium rosenbergii", "Decapoda", "Palaemonidae", "Crustacea", "Giant Freshwater Prawn", "Golda Chingri", "GBM_Freshwater,Sundarbans_Estuarine", "Core_Commercial", "LC", []),
    ("Macrobrachium malcolmsonii", "Decapoda", "Palaemonidae", "Crustacea", "Monsoon River Prawn", "Chatka Chingri", "GBM_Freshwater,Sundarbans_Estuarine", "Core_Commercial", "LC", []),
    ("Macrobrachium rude", "Decapoda", "Palaemonidae", "Crustacea", "Hairy River Prawn", "Goda Chingri", "GBM_Freshwater,Sundarbans_Estuarine", "Core_Commercial", "LC", []),
    ("Macrobrachium villosimanus", "Decapoda", "Palaemonidae", "Crustacea", "Dimua River Prawn", "Dimua Ichha", "GBM_Freshwater,Sundarbans_Estuarine", "Artisanal_Fishery", "LC", []),
    ("Macrobrachium lamarrei", "Decapoda", "Palaemonidae", "Crustacea", "Kuncho River Prawn", "Gura Chingri", "GBM_Freshwater", "Core_Commercial", "LC", []),
    ("Macrobrachium mirabile", "Decapoda", "Palaemonidae", "Crustacea", "Shortleg River Prawn", "Lotia Ichha", "GBM_Freshwater,Sundarbans_Estuarine", "Artisanal_Fishery", "LC", []),
    ("Exopalaemon styliferus", "Decapoda", "Palaemonidae", "Crustacea", "Roshna Prawn", "Roshna Chingri", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", ["Palaemon styliferus"]),
    ("Scylla olivacea", "Decapoda", "Portunidae", "Crustacea", "Orange Mud Crab", "Shila Kankra", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Scylla serrata", "Decapoda", "Portunidae", "Crustacea", "Giant Mud Crab", "Habba Kankra", "Sundarbans_Estuarine,Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Portunus pelagicus", "Decapoda", "Portunidae", "Crustacea", "Blue Swimming Crab", "Sataru Kankra", "Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Portunus sanguinolentus", "Decapoda", "Portunidae", "Crustacea", "Three-spot Swimming Crab", "Tin-fota Kankra", "Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Charybdis feriata", "Decapoda", "Portunidae", "Crustacea", "Crucifix Crab", "Kankra", "Bay_of_Bengal_Marine", "Artisanal_Fishery", "LC", []),
    ("Panulirus polyphagus", "Decapoda", "Palinuridae", "Crustacea", "Mud Spiny Lobster", "Chingri Lobster", "Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Panulirus ornatus", "Decapoda", "Palinuridae", "Crustacea", "Ornate Spiny Lobster", "Rongin Lobster", "Bay_of_Bengal_Marine", "Core_Commercial", "LC", []),
    ("Thenus orientalis", "Decapoda", "Scyllaridae", "Crustacea", "Flathead Lobster", "chapta Lobster", "Bay_of_Bengal_Marine", "Artisanal_Fishery", "LC", []),
)


def _deterministic_taxon_id(valid_name: str) -> int:
    """Create a deterministic positive synthetic taxon ID when external AphiaID/txid is unset."""
    digest = hashlib.sha256(valid_name.lower().encode("utf-8")).hexdigest()[:10]
    return 9_000_000 + (int(digest, 16) % 900_000)


def normalize_binomial(name: str) -> str:
    """Clean binomial string to 'Genus species' formatting."""
    parts = re.sub(r"[^A-Za-z\s-]", "", name.strip()).split()
    if len(parts) < 2:
        return ""
    genus = parts[0].capitalize()
    species = parts[1].lower()
    if species in ("sp", "spp", "cf", "aff"):
        return ""
    return f"{genus} {species}"


COMMERCIAL_CRUSTACEAN_FAMILIES: set[str] = {
    "Penaeidae",
    "Palaemonidae",
    "Portunidae",
    "Sergestidae",
    "Palinuridae",
    "Scyllaridae",
    "Solenoceridae",
    "Pandalidae",
}


def fetch_obis_bay_of_bengal_checklist(timeout: int = 25) -> list[dict[str, Any]]:
    """Fetch WoRMS-harmonized marine & brackish fishes + commercial decapods in the Bay of Bengal / Bangladesh EEZ.

    Uses OBIS v3 `/checklist` endpoint (`size=2000`, single JSON response per clade).
    Bounding polygon covers Bangladesh coastal waters, Sundarbans, and Northern Bay of Bengal (88.0E-92.8E, 18.0N-23.0N).
    """
    wkt = "POLYGON((88.0 18.0,92.8 18.0,92.8 23.0,88.0 23.0,88.0 18.0))"
    # WoRMS AphiaIDs: Actinopteri (843664), Elasmobranchii (10193), Decapoda (1130)
    groups = [
        (843664, "Actinopterygii", 2),
        (10193, "Elasmobranchii", 2),
        (1130, "Crustacea", 5),
    ]
    discovered: list[dict[str, Any]] = []
    session = requests.Session()
    session.headers.update({"User-Agent": "BengalFishRef/1.0 (FBD-NIB Bangladesh; scientific checklist reconciliation)"})

    for taxon_id, group_label, gcode in groups:
        url = "https://api.obis.org/v3/checklist"
        params = {"geometry": wkt, "taxonid": taxon_id, "size": 2000}
        try:
            resp = session.get(url, params=params, timeout=timeout)
            resp.raise_for_status()
            payload = resp.json()
            for item in payload.get("results", []):
                if (item.get("taxonRank") or "").lower() != "species":
                    continue
                fam = item.get("family") or "Incertae_sedis"
                if group_label == "Crustacea" and fam not in COMMERCIAL_CRUSTACEAN_FAMILIES:
                    continue
                raw_sci = item.get("acceptedNameUsage") or item.get("scientificName") or ""
                vname = normalize_binomial(raw_sci)
                if not vname:
                    continue
                syns: list[str] = []
                orig_sci = normalize_binomial(item.get("scientificName") or "")
                if orig_sci and orig_sci != vname:
                    syns.append(orig_sci)
                zones: list[str] = []
                if item.get("is_freshwater"):
                    zones.append("GBM_Freshwater")
                if item.get("is_brackish"):
                    zones.append("Sundarbans_Estuarine")
                if item.get("is_marine") or not zones:
                    zones.append("Bay_of_Bengal_Marine")
                iucn = item.get("redlist_category") or "NE"
                if iucn not in ("CR", "EN", "VU", "NT", "LC", "DD", "NE"):
                    iucn = "NE"
                tid = int(item.get("ncbi_id") or item.get("taxonID") or _deterministic_taxon_id(vname))
                discovered.append(
                    {
                        "taxon_id": tid,
                        "valid_name": vname,
                        "authority": item.get("scientificNameAuthorship") or "",
                        "common_name_en": "",
                        "common_name_bn": "",
                        "tax_group": group_label,
                        "order_name": item.get("order") or "Incertae_sedis",
                        "family_name": fam,
                        "genus_name": vname.split()[0],
                        "synonyms": syns,
                        "habitat_zones": ",".join(zones),
                        "aquaculture_tier": "Artisanal_Fishery",
                        "iucn_bd_status": iucn,
                        "genetic_code": gcode,
                    }
                )
            time.sleep(0.35)
        except Exception:
            continue
    return discovered


def fetch_inaturalist_bangladesh_taxa(timeout: int = 25) -> list[dict[str, Any]]:
    """Fetch vouchered inland freshwater & coastal fish species observed in Bangladesh (place_id=7154)."""
    session = requests.Session()
    session.headers.update({"User-Agent": "BengalFishRef/1.0 (FBD-NIB Bangladesh; scientific checklist reconciliation)"})
    # Bangladesh place_id=7154; Actinopterygii=47178, Elasmobranchii=47273
    queries = [
        (47178, "Actinopterygii", 2, "GBM_Freshwater"),
        (47273, "Elasmobranchii", 2, "Bay_of_Bengal_Marine"),
    ]
    discovered: list[dict[str, Any]] = []

    for inat_tid, group_label, gcode, default_zone in queries:
        for page in (1, 2):
            try:
                resp = session.get(
                    "https://api.inaturalist.org/v1/observations/species_counts",
                    params={
                        "place_id": 7154,
                        "taxon_id": inat_tid,
                        "hrank": "species",
                        "lrank": "species",
                        "per_page": 500,
                        "page": page,
                    },
                    timeout=timeout,
                )
                resp.raise_for_status()
                data = resp.json()
                results = data.get("results", [])
                if not results:
                    break
                for entry in results:
                    tx = entry.get("taxon", {})
                    vname = normalize_binomial(tx.get("name", ""))
                    if not vname:
                        continue
                    cs = (tx.get("conservation_status") or {}).get("status_name", "")
                    iucn_map = {
                        "critically endangered": "CR",
                        "endangered": "EN",
                        "vulnerable": "VU",
                        "near threatened": "NT",
                        "least concern": "LC",
                        "data deficient": "DD",
                    }
                    iucn = iucn_map.get(cs.lower(), "NE")
                    discovered.append(
                        {
                            "taxon_id": _deterministic_taxon_id(vname),
                            "valid_name": vname,
                            "authority": "",
                            "common_name_en": tx.get("preferred_common_name") or "",
                            "common_name_bn": "",
                            "tax_group": group_label,
                            "order_name": "Teleostei" if group_label == "Actinopterygii" else group_label,
                            "family_name": tx.get("iconic_taxon_name") or group_label,
                            "genus_name": vname.split()[0],
                            "synonyms": [],
                            "habitat_zones": default_zone,
                            "aquaculture_tier": "Artisanal_Fishery",
                            "iucn_bd_status": iucn,
                            "genetic_code": gcode,
                        }
                    )
                if len(results) < 500:
                    break
                time.sleep(0.35)
            except Exception:
                break
    return discovered


def build_master_checklist(
    conn: sqlite3.Connection,
    include_remote_reconciliation: bool = True,
    export_tsv: Path | None = None,
) -> dict[str, Any]:
    """Build and store the harmonized Bangladesh & Bay of Bengal Fish + Commercial Crustacean checklist."""
    merged: dict[str, dict[str, Any]] = {}

    if include_remote_reconciliation:
        for row in fetch_obis_bay_of_bengal_checklist():
            merged[row["valid_name"]] = row
        for row in fetch_inaturalist_bangladesh_taxa():
            v = row["valid_name"]
            if v in merged:
                existing_zones = set(merged[v]["habitat_zones"].split(","))
                new_zones = set(row["habitat_zones"].split(","))
                merged[v]["habitat_zones"] = ",".join(sorted(existing_zones | new_zones))
                if not merged[v]["common_name_en"] and row["common_name_en"]:
                    merged[v]["common_name_en"] = row["common_name_en"]
                if merged[v]["iucn_bd_status"] == "NE" and row["iucn_bd_status"] != "NE":
                    merged[v]["iucn_bd_status"] = row["iucn_bd_status"]
            else:
                merged[v] = row

    # Overlay authoritative curated FBD-NIB anchor metadata (highest precedence for BN names, tiers, synonyms)
    for (
        valid_name,
        order_name,
        family_name,
        tax_group,
        common_en,
        common_bn,
        zones,
        tier,
        iucn_bd,
        synonyms,
    ) in CURATED_BANGLADESH_ANCHOR_TAXA:
        v = normalize_binomial(valid_name)
        gcode = 5 if tax_group == "Crustacea" else 2
        existing = merged.get(v, {})
        combined_zones = set(zones.split(","))
        if existing.get("habitat_zones"):
            combined_zones |= set(existing["habitat_zones"].split(","))
        merged[v] = {
            "taxon_id": existing.get("taxon_id") or _deterministic_taxon_id(v),
            "valid_name": v,
            "authority": existing.get("authority", ""),
            "common_name_en": common_en or existing.get("common_name_en", ""),
            "common_name_bn": common_bn,
            "tax_group": tax_group,
            "order_name": order_name,
            "family_name": family_name,
            "genus_name": v.split()[0],
            "synonyms": synonyms,
            "habitat_zones": ",".join(sorted(z for z in combined_zones if z)),
            "aquaculture_tier": tier,
            "iucn_bd_status": iucn_bd,
            "genetic_code": gcode,
        }

    # Remove any synonym duplicates that appear under an older genus name in remote checklists
    all_synonyms_to_valid: dict[str, str] = {}
    for v, item in list(merged.items()):
        for syn in item.get("synonyms", []):
            ns = normalize_binomial(syn)
            if ns and ns != v:
                all_synonyms_to_valid[ns] = v
    for syn_name in all_synonyms_to_valid:
        if syn_name in merged:
            del merged[syn_name]

    ordered_rows = [merged[k] for k in sorted(merged.keys())]
    n_upserted = upsert_taxa(conn, ordered_rows)

    if export_tsv is not None:
        export_tsv.parent.mkdir(parents=True, exist_ok=True)
        with export_tsv.open("w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(
                [
                    "taxon_id",
                    "valid_name",
                    "tax_group",
                    "order_name",
                    "family_name",
                    "genus_name",
                    "common_name_en",
                    "common_name_bn",
                    "habitat_zones",
                    "aquaculture_tier",
                    "iucn_bd_status",
                    "genetic_code",
                    "synonyms",
                ]
            )
            for r in ordered_rows:
                writer.writerow(
                    [
                        r["taxon_id"],
                        r["valid_name"],
                        r["tax_group"],
                        r["order_name"],
                        r["family_name"],
                        r["genus_name"],
                        r["common_name_en"],
                        r["common_name_bn"],
                        r["habitat_zones"],
                        r["aquaculture_tier"],
                        r["iucn_bd_status"],
                        r["genetic_code"],
                        ";".join(r["synonyms"]),
                    ]
                )

    by_group: dict[str, int] = {}
    for r in ordered_rows:
        by_group[r["tax_group"]] = by_group.get(r["tax_group"], 0) + 1

    return {
        "total_species": n_upserted,
        "by_tax_group": by_group,
        "synonyms_mapped": len(all_synonyms_to_valid),
    }


def resolve_species_name_Map(conn: sqlite3.Connection) -> dict[str, dict[str, Any]]:
    """Build lookup dictionary mapping valid names AND historical synonyms -> taxon row dict."""
    lookup: dict[str, dict[str, Any]] = {}
    cur = conn.execute("SELECT * FROM taxa")
    for row in cur.fetchall():
        d = dict(row)
        lookup[d["valid_name"].lower()] = d
        for syn in json.loads(d["synonyms_json"] or "[]"):
            ns = normalize_binomial(syn).lower()
            if ns:
                lookup[ns] = d
    return lookup


def generate_dark_taxa_gap_report(conn: sqlite3.Connection) -> list[dict[str, Any]]:
    """Generate the per-species marker coverage and BDTPI prioritization matrix."""
    taxa_rows = [dict(r) for r in conn.execute("SELECT * FROM taxa ORDER BY valid_name").fetchall()]
    counts_cur = conn.execute(
        """
        SELECT resolved_taxon_id, marker, quality_tier, COUNT(*) AS n
        FROM sequences
        WHERE resolved_taxon_id IS NOT NULL
        GROUP BY resolved_taxon_id, marker, quality_tier
        """
    )
    gold_map: dict[int, dict[str, int]] = {}
    usable_map: dict[int, dict[str, int]] = {}
    for r in counts_cur.fetchall():
        tid = int(r["resolved_taxon_id"])
        m = str(r["marker"])
        tier = str(r["quality_tier"])
        n = int(r["n"])
        if tier == "Gold":
            gold_map.setdefault(tid, {})[m] = gold_map.get(tid, {}).get(m, 0) + n
        if tier in ("Gold", "Silver"):
            usable_map.setdefault(tid, {})[m] = usable_map.get(tid, {}).get(m, 0) + n

    report: list[dict[str, Any]] = []
    for tx in taxa_rows:
        tid = int(tx["taxon_id"])
        gm = gold_map.get(tid, {})
        um = usable_map.get(tid, {})
        bdtpi = compute_bdtpi(tx["iucn_bd_status"], tx["aquaculture_tier"], gm)
        is_dark_all = sum(um.values()) == 0
        is_dark_edna = (um.get("12S", 0) + um.get("16S", 0) + um.get("Mitogenome", 0)) == 0
        entry = {
            "taxon_id": tid,
            "valid_name": tx["valid_name"],
            "common_name_en": tx["common_name_en"],
            "common_name_bn": tx["common_name_bn"],
            "tax_group": tx["tax_group"],
            "family_name": tx["family_name"],
            "habitat_zones": tx["habitat_zones"],
            "aquaculture_tier": tx["aquaculture_tier"],
            "iucn_bd_status": tx["iucn_bd_status"],
            "bdtpi_score": bdtpi,
            "is_completely_dark": is_dark_all,
            "is_edna_12s_16s_dark": is_dark_edna,
        }
        for m in CORE_MARKERS:
            entry[f"{m}_gold"] = gm.get(m, 0)
            entry[f"{m}_usable"] = um.get(m, 0)
        report.append(entry)

    report.sort(key=lambda x: x["bdtpi_score"], reverse=True)
    return report
