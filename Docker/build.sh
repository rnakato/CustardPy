reponame=custardpy

tag=1.0.0
docker build -f Dockerfile.$tag -t rnakato/$reponame:$tag . #--no-cache
#docker save -o image_name.tar rnakato/$reponame:$tag
#singularity build -F $reponame.$tag.sif docker-archive://image_name.tar
docker push rnakato/$reponame:$tag
docker tag rnakato/$reponame:$tag rnakato/$reponame:latest
docker push rnakato/$reponame:latest
