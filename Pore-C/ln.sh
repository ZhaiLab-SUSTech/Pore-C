for i in `ls -d /public/home/lizw/task/pore_c/tools/*/`;do
mkdir ${i:36}
ln ${i}/* ${i:36}/
done

for a in `ls /public/home/lizw/task/pore_c/tools`;do
ln /public/home/lizw/task/pore_c/tools/${a} .
done
