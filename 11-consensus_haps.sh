cat ids.txt | while read line; do sed -i "s/>/>${line}_/" ${line}.20g085100.fa; done
