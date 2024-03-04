# Docker info for Week 4

command to build is:
```
docker1 build -f Dockerfile -t scrna2024 `pwd`
```

If you add commands to the end of the Dockerfile you can call 'build', and it will add a layer to previously built image, without having to redo all the previous steps.

Once you are happy with an image, save it to disk with `docker1 save -o image.tar imageID`
