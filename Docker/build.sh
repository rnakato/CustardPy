reponame=custardpy

tag=3.2.1
docker build -f Dockerfile.$tag -t rnakato/$reponame:$tag . #--no-cache
#apptainer build -F /work3/SingularityImages/$reponame.$tag.sif docker-daemon://rnakato/$reponame:$tag
#exit
docker push rnakato/$reponame:$tag
docker tag rnakato/$reponame:$tag rnakato/$reponame:latest
docker push rnakato/$reponame:latest
