#!/bin/bash

# Controlla se è stato passato un percorso come argomento
if [ -z "$1" ]; then
  echo "Usa: $0 <percorso>"
  exit 1
fi

# Memorizza il percorso passato come argomento
DIR=$1

# Controlla se il percorso esiste ed è una directory
if [ ! -d "$DIR" ]; then
  echo "Errore: $DIR non è una directory valida"
  exit 1
fi

# Modifica i permessi delle cartelle a 777
find "$DIR" -type d -exec chmod 777 {} \;

# Modifica i permessi dei file con estensione .x e .sh a 777
find "$DIR" -type f \( -name "*.x" -o -name "*.sh" \) -exec chmod 777 {} \;

# Modifica i permessi di tutti gli altri file a 666
find "$DIR" -type f ! \( -name "*.x" -o -name "*.sh" \) -exec chmod 666 {} \;

echo "Permessi modificati con successo nel percorso: $DIR"
