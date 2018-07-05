web: sh -c 'cd ./genviz/ && exec gunicorn genviz.wsgi --log-file -'
release: python genviz/manage.py migrate
