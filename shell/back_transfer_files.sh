rsync -av -e ssh --exclude-from 'shell/exclude_list.txt' lbenz@login.rc.fas.harvard.edu:~/Levis_Phd_Paper1-ExtraSims/fully_flexible/* ./fully_flexible
rm -r fully_flexible/logs
rm -r fully_flexible/model_objects