# 🚀 MAKE QUALITY GREAT AGAIN

**Module open source d'autoqualification de MNT (GEMAUT) — délivre une carte de précision associée au MNT (GEMAUT) **

Développé au Service de l'Imagerie Spatiale de l'IGN.

## 📖 Description

MAKE QUALITY GREAT AGAIN (MQGA) calcule un **masque de qualité / carte de précision** à partir :

- d'un **MNT** (DTM), typiquement produit par GEMAUT 
- d'un **MNS** (DSM)

Ce dépôt est le volet d'auto-qualification de [GEMAUT-pipeline](https://github.com/IGNF/GEMAUT-pipeline).

L'algorithme est décrit en détail dans cette article publié dans la *Revue Française de Photogrammétrie et de Télédétection (RFPT)* :

> 👉 [Lire l'article sur le site de la RFPT](https://rfpt.sfpt.fr/index.php/RFPT/article/view/739)

---

## 🏗️ Installation

### Méthode recommandée : environnement conda

```bash
# Cloner le dépôt
git clone https://github.com/IGNF/make-quality-great-again.git
cd make-quality-great-again

# Créer et activer l'environnement
conda env create -f mqga_env.yml
conda activate mqga_env

# Vérifier l'installation
python3 make_quality_great_again.py --help
```

---

## 🎯 Utilisation

### 1. Activer l'environnement

```bash
conda activate mqga_env
```

### 2. Aide du script

```bash
python3 make_quality_great_again.py --help
```

### 3. Exemples d'utilisation

#### **Cas standard (recommandé) avec interpolation hybride**

```bash
python3 make_quality_great_again.py \
    --mns /chemin/vers/MNS.tif \
    --mnt /chemin/vers/MNT.tif \
    --out /chemin/vers/masque_qualite.tif \
    --RepTra /chemin/vers/tmp_mqga \
    --cpu 8 \
    --interp hybrid \
    --clean
```

#### **Mode verbeux (DEBUG console)**

```bash
python3 make_quality_great_again.py \
    --mns /chemin/vers/MNS.tif \
    --mnt /chemin/vers/MNT.tif \
    --out /chemin/vers/masque_qualite.tif \
    --RepTra /chemin/vers/tmp_mqga \
    --cpu 8 \
    --interp hybrid \
    --verbose
```

#### **Avec détection de décrochage MNT**

```bash
python3 make_quality_great_again.py \
    --mns /chemin/vers/MNS.tif \
    --mnt /chemin/vers/MNT.tif \
    --out /chemin/vers/masque_qualite.tif \
    --RepTra /chemin/vers/tmp_mqga \
    --cpu 8 \
    --interp hybrid \
    --decrochage \
    --clean
```

Les zones où le MNT semble décrocher (MNS ≫ MNT, après seuil Z / morpho / volume / STD)
sont mises en **NoData** dans la carte de précision, et exportées en shapefile
`masque_qualite_zones_decrochage.shp` à côté de `--out`.

---

## 🔧 Paramètres

### Paramètres obligatoires

- `--mns` : MNS d'entrée (DSM)
- `--mnt` : MNT d'entrée (DTM)
- `--out` : Masque de qualité en sortie
- `--cpu` : Nombre de CPUs à utiliser
- `--RepTra` : Répertoire de travail temporaire

### Paramètres optionnels

- `--no` : Valeur NoData (défaut: `-9999`)
- `--per` : Percentile local (défaut: `0.05`)
- `--demiwinl` : Demie-taille en ligne de la fenêtre d'analyse (défaut: `50`)
- `--demiwinc` : Demie-taille en colonne de la fenêtre d'analyse (défaut: `50`)
- `--tile` : Taille de la tuile (défaut: `500`)
- `--pad` : Recouvrement entre tuiles (défaut: `50`)
- `--winavg` : Fenêtre de moyenne glissante finale (défaut: `50`)
- `--interp` : Méthode d'interpolation des NoData  
  (`hybrid` par défaut, ou `idw`)
- `--clean` : Vider le répertoire temporaire s'il existe déjà
- `--verbose` : Activer le niveau DEBUG sur la console
- `--decrochage` : activer la détection des zones de décrochage MNT
- `--seuilZ` : seuil Z positif MNS−MNT (m), défaut `10`
- `--seuilV` : seuil de volume `COUNT*MEAN`, défaut `100000`
- `--seuilSTD` : seuil d'écart-type (m), défaut `3`
- `--morph-radius` : rayon (px) de l'opening morphologique, défaut `5`

---

## 📁 Structure des données

### Entrées

- **MNS** et **MNT** au format GeoTIFF
- Ils peuvent avoir des **résolutions différentes** ; la carte de précision en sortie est calée sur la **grille du MNS**
- Si un pixel MNS ou MNT est NoData → pixel NoData dans la différence

### Sorties

- `--out` : masque de qualité / carte de précision LE90 estimée locale (GeoTIFF, résolution MNS)
- Avec `--decrochage` : shapefile `*_zones_decrochage.shp` (zones non garanties) ;
  ces pixels restent **NoData** dans `--out` (non interpolés)
- Fichier log : même chemin que `--out`, extension `.log`
- Fichiers temporaires dans `--RepTra` (dont `diff_MNS_MNT.tif`, `mask_decrochage.tif`)

### Organisation du traitement

1. Calcul de la différence MNS − MNT (`compute_mns_mnt_diff`)
2. (Optionnel) Détection décrochage → masque + shapefile ; masquage de la diff
3. Lecture des métadonnées (`GetInfo` → `RasterInfo`)
4. Calcul du nombre de tuiles (`CalculNombreDallesXY` → `TileGrid`)
5. Découpage en tuiles (`MakeDecoupage`)
6. Calcul de la carte de précision par tuile en parallèle (`DoParallel`)
7. Assemblage final (`Make_Assemblage_FINAL`)
8. Interpolation des NoData (`interpolate_nodata_*`, hors zones de décrochage)
9. Lissage final (`apply_moving_average`)

---

## 🩹 Interpolation hybride (`--interp hybrid`) — recommandée

Équivalent Python de :

```text
xingng -FB:2:C:50,1:1:50:1 -EM=-9999
```

### Principe

1. Identification des **trous connexes** (zones NoData)
2. Détection des **pixels de bord** (connexité 4)
3. Pour chaque trou :
   - calcul d'une constante `V_calc` : après exclusion des 50 % plus petites valeurs de bord, `V_calc` correspond au minimum des valeurs restantes
   - interpolation **IDW** locale sur les pixels de bord (rayon 50, poids 1)
   - combinaison selon la distance au bord : proche → IDW, centre → `V_calc`

Cette méthode s'inspire de la méthode `xingng`: elle remplit bien les grands trous tout en restant raisonnablement dans son temps d'exécution.

---

## 🗂️ Structure du dépôt

```text
MAKE_QUALITY_GREAT_AGAIN/
├── make_quality_great_again.py   # Point d'entrée CLI (wrapper)
├── mqga/
│   ├── cli.py                    # Arguments, logging, orchestration
│   ├── io_raster.py              # I/O raster, diff MNS−MNT
│   ├── decrochage.py             # Détection zones MNS ≫ MNT
│   ├── quality_mask.py           # Percentile local / masque
│   ├── tiling.py                 # Découpage et parallèle
│   ├── mosaic.py                 # Assemblage des tuiles
│   └── interpolate.py            # Bouchage NoData + lissage
├── mqga_env.yml                  # Environnement conda
└── README.md
```

---

## Contribution

1. Fork le projet
2. Créer une branche feature (`git checkout -b feature/AmazingFeature`)
3. Commit les changements (`git commit -m 'Add AmazingFeature'`)
4. Push vers la branche (`git push origin feature/AmazingFeature`)
5. Ouvrir une Pull Request

---

## Licence

Ce projet est sous la licence [LICENSE](LICENSE).
