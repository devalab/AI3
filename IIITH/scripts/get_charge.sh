set -x
export PYTHONUNBUFFERED=1
<<<"addcharge" chimera --nogui --nostatus $1 > charge.dat & pid=$!
sleep 2
kill -INT "$pid"
pkill -f "antechamber"