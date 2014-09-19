#!/usr/bin/env python
from PyQt4.Qt import QApplication, QListView, QStandardItemModel, QStandardItem
import sys
from obspy import readEvents

cat = readEvents("/home/megies/svn/mixed/katalog_unterhaching_stand_2013-04-30.intern.xml")


# Create a Qt application
app = QApplication(sys.argv)

# Our main window will be a QListView
list = QListView()
list.setWindowTitle('Example List')
list.setMinimumSize(600, 400)

# Create an empty model for the list's data
model = QStandardItemModel(list)

## # Add some textual items
## foods = [
##     'Cookie dough', # Must be store-bought
##     'Hummus', # Must be homemade
##     'Spaghetti', # Must be saucy
##     'Dal makhani', # Must be spicy
##     'Chocolate whipped cream' # Must be plentiful
## ]
## 
## for food in foods:
##     # create an item with a caption
##     item = QStandardItem(food)
## 
##     # add a checkbox to it
##     item.setCheckable(True)
## 
##     # Add the item to the model
##     model.appendRow(item)

for pick in cat[-40].picks:
     # create an item with a caption
     item = QStandardItem(str(pick.resource_id))
 
     # add a checkbox to it
     item.setCheckable(True)
     item.setEditable(False)
     item.setTristate(True)
 
     # Add the item to the model
     model.appendRow(item)

# Apply the model to the list view
list.setModel(model)

# Show the window and run the app
list.show()
app.exec_()
