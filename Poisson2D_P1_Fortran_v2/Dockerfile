# Environnement Python sur une base Debian.
# L'étiquette fixe la version mineure de Python et la famille du système.
FROM python:3.11-slim-bookworm

# Limiter les bibliothèques concernées à un seul thread.
ENV OMP_NUM_THREADS=1
ENV OPENBLAS_NUM_THREADS=1

# Installer le compilateur Fortran et Make.
# Le nettoyage réduit la place occupée par les listes de paquets.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       gfortran \
       make \
    && rm -rf /var/lib/apt/lists/*

# Créer un utilisateur non privilégié pour compiler et exécuter le projet.
RUN groupadd --gid 1000 app \
    && useradd --uid 1000 --gid app --create-home app

# Dossier principal du projet dans le conteneur.
WORKDIR /app

# Installer d'abord les dépendances Python.
# Si les sources changent mais pas requirements.txt,
# Docker peut réutiliser cette étape depuis son cache.
COPY requirements.txt .

RUN python -m pip install --no-cache-dir -r requirements.txt \
    && python -m pip check

# Autoriser l'utilisateur app à écrire dans le dossier du projet.
RUN chown app:app /app

# Copier les sources, en excluant les chemins de .dockerignore.
COPY --chown=app:app . .

# Les opérations suivantes ne sont plus exécutées en tant que root.
USER app

# Compiler et exécuter les tests, puis compiler le solveur optimisé.
# Si un test échoue, la construction de l'image échoue aussi.
# Ne pas conserver les résultats de tests dans l'image finale.
RUN make clean \
    && make test \
    && make \
    && rm -rf results \
    && mkdir -p results

# Commande par défaut : calculs du cas 1 et création des images.
# exec remplace le shell par le processus Python.
CMD exec python python/plot_results.py
