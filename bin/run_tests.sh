compilation(){
    # Exception : petsc clamped bar : compilation with mpiCCE
    if [ "$1" = "petsc_clamped_bar" ] && [ "$2" = "cpp=clang" ]; then
        echo $1 cannot compile with clang.
        return
    elif [ "$1" = "petsc_clamped_bar" ]; then
        echo $1 : Compilation with mpiCC and $2
        option="cpp=mpiCC $2"
    elif [ "$#" = 1 ]; then
        echo $1 : Compilation
    else
        echo $1 : Compilation with option \"$2\"
        option="$2"
    fi
    # Others can be compiled with gcc or clang.
    scons -j 8 $option --config=force > .test-compilation.log 2>&1

    if [ "$?" = "0" ]; then
        echo Compilation - Ok
    else
        echo -e "\033[1;31m  $1 Compilation - Fail \033[0m"
        cat .test-compilation.log
    fi
}

clean(){
    scons -c > /dev/null
}


run(){
    # Exception : the forward must be run first with different configuration.
    if [ -f "./forward" ]; then
        echo $1 : forward : Run
        ./forward configuration/truth.lua > .forward-run.log 2>&1
            if [ "$?" = "0" ]; then
                echo Run - Ok
            else
                echo -e '\033[1;31m  Run - Fail \033[0m'
                cat .forward-run.log
            fi
    fi
    # Exception : the generate_observation must be run second with different configuration.
    if [ -f "./generate_observation" ]; then
        echo $1 : generate_observation : Run
        ./generate_observation configuration/truth.lua > .generate_observation-run.log 2>&1
            if [ "$?" = "0" ]; then
                echo Run - Ok
            else
                echo -e '\033[1;31m  Run - Fail \033[0m'
                cat .generate_observation-run.log
            fi
    fi
    # All other executables can be run with the correct configuration.
    for execable in *; do
        if [ "$execable" = "forward" ] || [ "$execable" = "generate_observation" ]; then
            continue
        fi
        if [ -x "$execable" ] && [ ! -d "$execable" ]; then
            echo $1 : $execable : Run
            if [ "$execable" = "run" ];then
                config=""
            else
                config="configuration/assimilation.lua"
            fi
            ./$execable $config > .$execable-run.log 2>&1
            if [ "$?" = "0" ]; then
                echo Run - Ok
            else
                echo -e '\033[1;31m  Run : $execable : Fail \033[0m'
                cat .$execable-run.log
            fi
        fi
    done
}



if [ "$1" = "help" ]; then
echo Usage :
echo \"$0\" will display only failed instruction, everything will be logged. Tests will be compiled with gcc.
echo \"$0 all\" tests will be compiled with gcc, clang and C++0x.
echo \"$0 examples\" examples will be compiled and run.
echo \"$0 norun\" tests or examples will be compiled but not run.
exit
fi


options[0]=" "
if [ "$1" = "all" ] || [ "$2" = "all" ] || [ "$3" = "all" ]; then
    options[1]="std=2011"
    options[2]='cpp=clang'
fi


# Reading command line options to know if examples must be compiled.
if [ "$1" = "examples" ] || [ "$2" = "examples" ] || [ "$3" = "examples" ]; then
    cd ../example/
    echo Examples tests
else
    cd ../test/unit/
    echo Unit tests
fi

# Loops on directories (inside /test/unit or /example).
for dir in *; do
    cd $dir

    # Loops on options (e.g : C++0x, Clang...).
    for i in ${!options[@]}; do
        clean
        compilation $dir ${options[i]}
        if [ "$1" != "norun" ] && [ "$2" != "norun" ] && [ "$3" != "norun" ] && [ "$dir" != "template" ]; then
            run $dir
        fi
    done
    echo " "
    cd ..
done 
cd ..

