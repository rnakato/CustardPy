reponame=custardpy

tag=2.0.0
docker build -f Dockerfile.$tag -t rnakato/$reponame:$tag . #--no-cache
#echo "docker save -o $reponame-$tag.tar rnakato/$reponame:$tag"
#docker save -o $reponame-$tag.tar rnakato/$reponame:$tag
#singularity build -F /work3/SingularityImages/$reponame.$tag.sif docker-archive://$reponame-$tag.tar
#exit

docker push rnakato/$reponame:$tag
docker tag rnakato/$reponame:$tag rnakato/$reponame:latest
docker push rnakato/$reponame:latest
