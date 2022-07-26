for tag in 0.1.0
do
    docker build -f Dockerfile.$tag -t rnakato/custardpy:$tag .
    docker push rnakato/custardpy:$tag
    docker build -f Dockerfile.$tag -t rnakato/custardpy .
    docker push rnakato/custardpy
done
