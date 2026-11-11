foo=(a1 b2 c4)
foo=("${foo[@]}" d7) # push element 'd'

foo+=(e3)

foo[${#foo[*]}]="f10" # push element 'f'

for i in "${foo[@]}"; do echo "$i" ; done

echo ${foo[@]}

echo now testing shell script list input into python

python3 arg_test.py ${foo[@]}