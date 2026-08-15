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

Les zones où le MNT semble décrocher sont mises en **NoData** dans la carte de précision et exportées en shapefile
`masque_qualite_zones_decrochage.shp` à côté de `--out`.

---

## 🔧 Paramètres

Alignés sur `python3 make_quality_great_again.py --help`.

### Options obligatoires

- `--mnt` : chemin vers le MNT (DTM)
- `--mns` : chemin vers le MNS (DSM)
- `--out` : masque de qualité / carte de précision en sortie
- `--RepTra` : répertoire de travail (fichiers temporaires)
- `--cpu` : nombre de CPU disponibles

### Options générales

- `--reso` : résolution de travail en mètres (ex. `1` ou `4`) ; si absent = résolution du MNS
- `--no` : valeur NoData (défaut: `-9999`)
- `--clean` : vider le répertoire temporaire s'il existe déjà
- `--verbose` : activer le niveau DEBUG sur la console

### Découpage en tuiles

- `--tile` : taille d'une tuile en pixels (défaut: `500`)
- `--pad` : recouvrement entre tuiles en pixels (défaut: `50` ; doit être `>= --demiwin`)

### Calcul de ε_stat (basée sur la partie négative de l'histogramme des différences MNS−MNT Cf. Papier RFPT)

- `--bias` : biais `|b|` ajouté à ε_stat : `ε_stat ← |b| + ε_stat` (défaut: `0`)
- `--stat` : ε_stat calculée par `percentile` (défaut) ou `mad` sur une fenêtre locale 
- `--per` : défaut: `0.10` ; ignoré si `mad`
- `--mad-k` : facteur `k` dans `ε_stat = |b| + k·MAD` si `--stat mad` (défaut: `2.44` ≈ LE90 ; ignoré si `percentile`)
- `--demiwin` : demi-taille (pixels) de la fenêtre d'analyse (défaut: `50`)
- `--min-valid` : effectif minimal pour calculer la stat [STANAG] ; seuil effectif = `max(min-valid, min-valid-pct% de la fenêtre)` (défaut: `167`)
- `--min-valid-pct` : taux minimal (%) de pixels valides dans la fenêtre d'analyse; valide = à la fois négatif dans la différence MNS−MNT && `!=` NoData (défaut: `10`) ; `0` = désactive l'option

### Interpolation des NoData (trous)

- `--interp` : `hybrid` (défaut, recommandé) ou `idw`
  - **`hybrid`** : pour chaque pixel du trou, `ε₀` = transition IDW (près du bord) → P90(ε_bord) (au centre) ; puis rampe optionnelle (`--hole-alpha`)
  - **`idw`** : IDW seul (pas de P90 / pas de rampe) ; `--hole-alpha` et `--hole-lambda` sont ignorés
- `--hole-alpha` : ça pénalise avec une rampe plus on s'enfonce dans le trou ; pente de la rampe `α·d·Δ` (m d'incertitude par m de distance au bord) ; formule `ε = min(ε₀ + α·d·Δ, λ·P90)` ; `0` = sans rampe (`ε = ε₀`) ; défaut `0.01`
- `--hole-lambda` : plafond de la rampe `ε ≤ λ·P90(ε_bord)` si `--hole-alpha > 0` (ignoré sinon) ; défaut `1.5`

### Lissage final

- `--winavg` : taille (pixels) de la moyenne glissante finale (défaut: `50`)

### Zones de non garantie (décrochage MNT)

- `--decrochage` : détecter les zones où le MNT semble décrocher (MNS ≫ MNT) ; NoData dans la carte + shapefile à côté de `--out`
- `--seuilZ` : écart minimal MNS−MNT (m) pour qu'une une zone soit candidate (défaut: `10`)
- `--seuilV` : volume minimal d'une zone (nb pixels × écart moyen) pour la retenir (défaut: `100000`)
- `--seuilSTD` : écart-type minimal (m) de MNS−MNT dans la zone pour la retenir (défaut: `3`)
- `--morph-radius` : nettoyage morphologique (rayon en pixels) : enlève le bruit / petites taches (défaut: `5`)

---

## 📁 Structure des données

### Entrées

- **MNS** et **MNT** au format GeoTIFF
- Ils peuvent avoir des **résolutions différentes** ; par défaut la carte de précision a la résolution du MNS 
- Avec `--reso` : MNS et MNT sont rééchantillonnés à la résolution de l'utilisateur 

### Sorties

- `--out` : masque de qualité / carte de précision LE90 estimée localement 
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

### Sans / avec rampe

- **Sans rampe** (`--hole-alpha 0`) : on garde `ε₀` uniquement
- **Avec rampe** (défaut `α = 0.01`) :
  `ε = min(ε₀ + α·d·Δ, λ·P90)`  
  avec `λ` = `--hole-lambda` (défaut `1.5`) et `P90` = P90 des ε de bord du trou

### Principe:

Pour chaque trou NoData :

1. Identification du trou (composante connexe) et des pixels de bord valides 
2. Au centre : `V_calc = P90(ε_bord)`
3. Pour chaque pixel du trou, avec `d` = distance (px) au bord le plus proche et rayon fixe `50` :
   - `K = min(1, d / 50)`
   - si `d > 50` : on affecte `V_calc = P90(ε_bord)`
   - si `d ≤ 50` mélange avec un `V_IDW classique` : `ε₀ = (1 - K) · V_IDW + K · V_calc`

| Distance au bord | Comportement de `ε₀` |
|------------------|----------------------|
| `d ≈ 0` (bord) | surtout **IDW** |
| `0 < d < 50` | **mélange** IDW ↔ P90 |
| `d ≥ 50` (cœur du trou) | **P90** seul |

Puis, si `--hole-alpha > 0`, pénalité distance + plafond `λ·P90` (voir ci-dessus).

Inspirée de `xingng -FB:2:C:50,1:1:50:1 ...` (IDW bord + constante au centre). A l'usage, ça remplit bien les trous, tout en restant raisonnable en termes de temps d'exécution.

Ne touche pas les zones `--decrochage` qui restent en NoData.

---

## 🗂️ Structure du dépôt

```text
MAKE_QUALITY_GREAT_AGAIN/
├── make_quality_great_again.py   # Point d'entrée CLI (wrapper → mqga.cli)
├── mqga/
│   ├── __init__.py
│   ├── cli.py                    # Arguments, logging, orchestration
│   ├── io_raster.py              # I/O raster, diff MNS−MNT, --reso
│   ├── decrochage.py             # Zones de non-garantie (MNS ≫ MNT)
│   ├── quality_mask.py           # Stat locale (percentile / MAD) + min-valid
│   ├── tiling.py                 # Découpage en tuiles + parallèle
│   ├── mosaic.py                 # Assemblage des tuiles
│   └── interpolate.py            # Bouchage NoData (hybrid/idw) + lissage
├── mqga_env.yml                  # Environnement conda
├── LICENSE
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
