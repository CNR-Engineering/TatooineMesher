Mailleur Tatooine
=================

Interpolateur/mailleur à partir de profils en travers et de lignes de contraintes

Langage : Python 3

## Description

* **Interpolateur** linéaire à l'aide de lignes directrices 2D :
    * densificateur de profils 1D : ajout de profils intermédiaires interpolés
    * génération d'un semis de points 2D entre profils en travers
* **Mailleur 2D** : triangulation Delaunay contrainte avec comme :
    * sommets : ceux issus de l'interpolateur
    * arêtes contraintes : celles qui sont le long des profils en travers et des lignes de contraintes

## Pré-requis

###Fichier de profils en travers

Format        | PK  | Nom | Couches
-------------------------------------
ST            | Oui | Non | Non
i3s           |     | Non | Non
shp POINTZ    | ?   | ?   | Oui
shp POLYLINEZ | ?   | ?   | Non

/!\ Le PK doit être unique !
TODO: S'il n'est pas dans le fichier d'entrée, alors le calculer avec l'axe hydraulique

* Lignes (i3s, shp, ST)
    * les profils ne sont pas nécessairement ordonnées, c'est l'axe hydraulique qui permet de les ré-ordonner
    * tous les profils sont décrits dans le même sens (rive gauche à droite ou inversement)
* Points (shp)
    * tous les profils sont décrits dans le même sens (rive gauche à droite ou inversement)

### Fichier de lignes de contraintes
* formats possibles : i2s ou shp
* les lignes de contraintes ne se croissent pas
* les lignes sont toutes orientées dans le même sens que l'axe hydraulique

### Fichier axe hydraulique
* formats possibles : i2s ou shp
* une seule polyligne orientée de l'amont vers l'aval
* elle intersecte tous les profils (et les épis)
* (l'axe hydraulique n'est pas considéré comme une ligne de contrainte et peut donc intersecter des lignes de contrainte)

## Fonctionnalités

* intersection entre les profils et les lignes :
    * la recherche des intersections est élargie (on tolère une intersection fictive s'ils sont distants d'une distance inférieure à <code>args.dist_max</code>, si cette variable n'est pas ''None'')
    * ces intersections définissent des lits (entre les profils amont et aval et entre les limites gauche et droite)
* écriture des fichiers de maillage en : t3s (BlueKenue), LandXML ou slf (Telemac)

## TODO
* avoid pandas dependency
* speed-up execution
* add option to project per `Lit` for all profiles
* read lignes contraintes from ST files