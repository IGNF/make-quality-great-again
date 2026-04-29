# MAKE QUALITY GREAT AGAIN

**Module open source d'autoqualification de GEMAUT**

Ce dépôt contient un script Python (`make_quality_great_again.py`) qui
calcule un **masque de qualité** à partir d'une **différence DSM/DTM**

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

    usage: make_quality_great_again.py [-h] [-diff DIFF] [-out OUT] [-no NO] [-per PER] [-demiwinl DEMIWINL] [-demiwinc DEMIWINC] [-tile TILE] [-pad PAD] [-RepTra REPTRA] [-cpu CPU]
                                       [-winavg WINAVG] [-interp {griddata,idw,idw_old,window,linearnd,fast,hybrid}] [-clean]

    MAKE QUALITY GREAT AGAIN ALL ZONE - Version Spatialisée & Parallélisée

    options:
      -h, --help            show this help message and exit
      -diff DIFF            Différence DSM/DTM en entrée
      -out OUT              Masque de Qualité en sortie
      -no NO                Valeur de No Data
      -per PER              Valeur de percentile
      -demiwinl DEMIWINL    Demie-taille en ligne de la fenêtre d'analyse
      -demiwinc DEMIWINC    Demie-taille en colonne de la fenêtre d'analyse
      -tile TILE            Tile / Taille de la tuile
      -pad PAD              Pad / Recouvrement entre tuiles
      -RepTra REPTRA        Répertoire de Travail
      -cpu CPU              Nombre de CPU diponibles
      -winavg WINAVG        Taille de la fenêtre glissante pour la moyenne (par défaut 50x50)
      -interp {griddata,idw,idw_old,window,linearnd,fast,hybrid}
                            Méthode d'interpolation pour les pixels nodata (défaut: idw - optimisé et parallélisé, hybrid = méthode hybride xingng)
      -clean                Supprimer le contenu du répertoire temporaire s'il existe déjà

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

1.  Lecture des métadonnées de l'image (`GetInfo`)
2.  Calcul du nombre de tuiles (`CalculNombreDallesXY`)
3.  Découpage en tuiles (`MakeDecoupage`)
4.  Calcul du masque par tuile en parallèle (`DoParallel`)
5.  Assemblage final des tuiles (`Make_Assemblage_FINAL`)
6.  Interpolation des NoData (`interpolate_nodata_*`)
7.  Lissage final (`apply_moving_average`)

------------------------------------------------------------------------

## Licence

Ce projet est sous licence [LICENSE](LICENSE).



