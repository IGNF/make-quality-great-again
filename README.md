# MAKE QUALITY GREAT AGAIN

**Module open source d'autoqualification de GEMAUT**

Ce dépôt contient un script Python (`make_quality_great_again.py`) qui
calcule un **masque de qualité** à partir d'un **MNS** et d'un **MNT**
(la différence MNS − MNT est calculée automatiquement).

------------------------------------------------------------------------

# Fonctionnalités principales

-   Découpage en **tuiles avec recouvrement** pour traiter de grands
    rasters
-   **Traitement parallèle (multiprocessing)** sur plusieurs CPU
-   Calcul d'un **masque de qualité par fenêtre glissante** (percentile
    sur la distribution locale)
-   **Assemblage automatique** des tuiles avec gestion du recouvrement
-   **Interpolation des trous (NoData)** avec plusieurs méthodes :
    -   `griddata`
    -   `idw` (IDW vectorisé et parallélisé)
    -   `idw_old`
    -   `window`
    -   `linearnd`
    -   `fast`
    -   `hybrid` (équivalent à `xingng -FB:2:C:50,1:1:50:1`)
-   **Moyenne glissante finale** pour lisser le masque

------------------------------------------------------------------------

# Installation typique

Créer l'environnement conda et lancer l'outil :

``` bash
conda env create -f mqga_env.yml
conda activate mqga_env

python3 make_quality_great_again.py --help
```

Vous devriez obtenir ceci :

    usage: make_quality_great_again.py [-h] --mns MNS --mnt MNT -out OUT [-no NO] [-per PER]
                                       [-demiwinl DEMIWINL] [-demiwinc DEMIWINC] [-tile TILE]
                                       [-pad PAD] -RepTra REPTRA -cpu CPU [-winavg WINAVG]
                                       [-interp {griddata,idw,idw_old,window,linearnd,fast,hybrid}]
                                       [-clean]

    options:
      --mns / -mns MNS      Chemin vers le MNS (DSM)
      --mnt / -mnt MNT      Chemin vers le MNT (DTM)
      -out OUT              Masque de Qualité en sortie
      -RepTra REPTRA        Répertoire de Travail
      -cpu CPU              Nombre de CPU disponibles
      ...

Exemple :

``` bash
python3 make_quality_great_again.py \
  --mns /chemin/mns.tif \
  --mnt /chemin/mnt.tif \
  -out /chemin/masque_qualite.tif \
  -RepTra /chemin/tmp_mqga \
  -cpu 8 \
  -interp hybrid \
  -clean
```

------------------------------------------------------------------------

✈️ **Et le voyage peut commencer !**

------------------------------------------------------------------------

# Méthode hybride (`-interp hybrid`) pour le bouchage de trous (recommnadée)

La méthode `hybrid` implémente en Python un équivalent de :
    xingng -FB:2:C:50,1:1:50:1 -EM=-9999
(Pour celles et ceux qui ont  la réf)

## Principe

1.  Identification des **trous connexes (zones NoData)**
2.  Détection des **pixels de bord du trou** (connexité 4)
3.  Pour chaque trou :

### Étape 1 --- Calcul d'une constante `V_calc`

-   récupération des pixels de bord
-   suppression des **50 % plus petites valeurs**
-   prise du **minimum des restantes**

### Étape 2 --- Interpolation locale

Interpolation **IDW sur les pixels de bord**

-   rayon : **50 pixels**
-   poids : **1**

### Étape 3 --- Combinaison des deux

Selon la **distance au bord du trou** :

  Position         Valeur utilisée
  ---------------- ----------------------
  proche du bord   interpolation locale
  centre du trou   constante `V_calc`

Cette méthode remplit mieux les grands trous tout en restant
**raisonnablement rapide**.

------------------------------------------------------------------------

# Organisation du traitement

1.  Calcul de la différence MNS − MNT sur la **grille du MNS**
    (`compute_mns_mnt_diff` ; le MNT est rééchantillonné si nécessaire)
2.  Lecture des métadonnées de l'image (`GetInfo`)
3.  Calcul du nombre de tuiles (`CalculNombreDallesXY`)
4.  Découpage en tuiles (`MakeDecoupage`)
5.  Calcul du masque par tuile en parallèle (`DoParallel`)
6.  Assemblage final des tuiles (`Make_Assemblage_FINAL`)
7.  Interpolation des NoData (`interpolate_nodata_*`)
8.  Lissage final (`apply_moving_average`)

------------------------------------------------------------------------

## Licence

Ce projet est sous licence [LICENSE](LICENSE).



