#########################################################################
# File Name: make.sh
# Author: Sun Wenliang
# mail: 11930830@mail.sustech.edu.cn
# Created Time: 2021年01月08日 星期五 15时52分39秒
#########################################################################
#!/bin/bash
cd ..
make cleanall
make
cd test
bash run_test.sh
