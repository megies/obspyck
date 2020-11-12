qt_designer: qt_designer.ui
	pyuic5 qt_designer.ui > obspyck/qt_designer.py
	# for whatever reason we have to manually fix the local import..
	sed --in-place 's#from util import QMplCanvas#from .util import QMplCanvas#' obspyck/qt_designer.py

# making a source distribution release (commented out because the git clean
# removes all local changes and untracked files)
#
# sdist:
# 	git checkout <tag-name>
# 	git clean -fdx
# 	umask 0022 && chmod -R a+rX . && python setup.py sdist --format=zip
# 	sha256sum dist/*
# 	twine upload dist/*
