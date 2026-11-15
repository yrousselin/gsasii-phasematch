@echo off
REM Installe la recherche de phases (search-match) dans GSAS-II (Windows)
REM Installs the phase search (search-match) plugin into GSAS-II (Windows)
set "DEST=%USERPROFILE%\.GSASII\exports"
if not exist "%DEST%" mkdir "%DEST%"
copy /Y "%~dp0G2export_PhaseMatch.py" "%DEST%\" >nul
copy /Y "%~dp0pm_*.py" "%DEST%\" >nul
copy /Y "%~dp0csd_pm_worker.py" "%DEST%\" >nul
echo.
echo Fichiers copies dans / Files copied to: %DEST%
echo Construisez ensuite l index de references (voir LISEZMOI.md / README.md).
echo Then build the reference index (see README.md / LISEZMOI.md).
echo Redemarrez GSAS-II / Restart GSAS-II:
echo   Export - Entire project as - (user) Phase search (search-match)
echo.
pause
