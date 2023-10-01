reponame=custardpy

tag=1.4.1
docker build -f Dockerfile.$tag -t rnakato/$reponame:$tag . #--no-cache

#docker save -o $reponame-$tag.tar rnakato/$reponame:$tag
#singularity build -F $reponame.$tag.sif docker-archive://$reponame-$tag.tar
#exit
docker push rnakato/$reponame:$tag
docker tag rnakato/$reponame:$tag rnakato/$reponame:latest
docker push rnakato/$reponame:latest
