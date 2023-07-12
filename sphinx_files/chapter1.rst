############
Sphinx メモ
############

******************
RTDスタイルの導入
******************
* pip で sphinx_rtd_theme を入れる。
* conf.pyで以下のように編集

.. code-block:: python

    import sphinx_rtd_theme
    html_theme = 'sphinx_rtd_theme'
    html_static_path = ['_static']
    html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]


***********
ロゴの導入
***********
* ロゴを描いてPNG形式で保存する。
* _static/ ディレクトリにロゴのpngファイルを入れる。今回はLOGO1.PNG
* conf.pyで以下を追記::

    html_logo = './_static/LOGO1.PNG'

* html_theme_options = {}の中をいじるともっと色々カスタマイズできるっぽい。

***********
toctree
***********
* rst拡張子を省いて並べて書くことで、指定したページへの参照が可能。
* 相対パスを使うこと。

****************************
ソースコードマニュアルの導入
****************************
* ディレクトリ構成は、doc/とsrc/とsphinx_files/が同じディレクトリ下にある状態
* sphinx_fileにこれまでのindex.rstなどがある
* src/ の下でコマンド実行時::

    sphinx-apidoc -F -H ProjectName -A AuthorName -V VersionNumber -o ../sphinx_files .

とすれば、既存のindex.rstなどは消されず、<module>.pyのマニュアルを呼ぶ<module>.rstがsphinx_files/に、<module>.htmlはdocs/に作られる。
* すでにあるconf.pyに対して

.. code-block:: python

    import os
    import sys
    sys.path.insert(0, '<Absolute PATH to src>')
    ~~~
    extensions = [
     'sphinx.ext.autodoc',
     'sphinx.ext.viewcode',
     'sphinx.ext.todo',
     'sphinx.ext.napoleon'
    ]

と編集してやれば良い。

