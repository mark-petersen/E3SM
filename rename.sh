for f in *.F90; do
    cp -- "$f" "${f%.F90}_single.F90"
done

