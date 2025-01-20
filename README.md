# macrotools
tools for polymer modeling and building


## Установка Pysimm
пререквизиты: lammps, conda (numpy, scipy, matplotlib)

Установка:

1. Копируем исходный код:

`git clone https://github.com/polysimtools/pysimm.git`

2. Создаем conda-окружение:

`conda create -n pysimm numpy scipy matplotlib`

`conda activate pysimm`

3. Устанавливаем исходный код:

`pip install -e /path/to/pysimm/folder`

Готово.

Если все было сделано правильно, теперь Пайсим всегда будет виден внутри окружения. Проверить можно командами:

`from pysimm import system, lmps, forcefield`

`from pysimm.apps.random_walk import random_walk, copolymer`

Перед запуском расчетов нужно добавить путь к исполняемому файлу Лампса в переменную `LAMMPS_EXEC`:

`export LAMMPS_EXEC=/path/to/lammps/executable`

## Построение ячейки аморфного полимера

..



