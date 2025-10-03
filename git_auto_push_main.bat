@echo off
cd /d %~dp0
echo Adding all changes to the Git repository...

git add .

echo Committing the changes...
set /p commitMsg="Enter commit message: "
git commit -m "%commitMsg%"

echo Pushing changes to remote...
git push origin main

echo Done!
pause
