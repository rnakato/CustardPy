#docker compose -f compose.yaml build r
docker compose -f compose.yaml build juicer
docker compose -f compose.yaml build custardpy

#exit
reponame=custardpy
tag=3.4.1
#apptainer build -F /work3/SingularityImages/$reponame.$tag.sif docker-daemon://rnakato/$reponame:$tag
#exit
docker push rnakato/$reponame:$tag
docker tag rnakato/$reponame:$tag rnakato/$reponame:latest
docker push rnakato/$reponame:latest
