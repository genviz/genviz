# Genviz
Herramienta de visualización de genes. Miniproyecto de desarrollo de software USB

## Instalación
Para instalar la herramienta en local se deben instalar las dependencias de Python y Javascript
```
# Crear y activar ambiente virtual
python3 -m venv env
source env/bin/activate

# Instalar paquetes de Python
pip install -r requirements.txt

# Instalar paquetes de JS
cd genviz
npm postinstall

# Correr migraciones
python genviz/manage.py migrate
```

## Deploy
La aplicación está siendo desplegada automáticamente a Heroku (https://usb-genviz.herokuapp.com) cada vez que se hace un commit a master.

## Generar diagrama de clases
Para generar de forma automatica un diagrama de clases a partir del modelo de Django, se debe hacer lo siguiente
```
sudo apt-get install graphviz
python genviz/manage.py graph_models app -o "docs/Class Diagram.png"
```