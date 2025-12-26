sed 's|void |void Propagation_param_t::|g' "$1" > "$2"
sed  's|int |int Propagation_param_t::|g' "$1" > "$2"
sed  's|double |real_t Propagation_param_t::|g' "$1" > "$2"
