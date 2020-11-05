cd "small/";

for d in */ ; do
    cd "$d";
    ./run.sh;
    cd "..";
done

cd "..";

for d in */ ; do
	if [[ "$d" != small* ]];
	then
    	cd "$d";
    	./run.sh;
    	cd "..";
    fi
done
