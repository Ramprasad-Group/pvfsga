# Test DB Seeds
In this folder is JSON that will be read into the test database. Each json file
corresponds to one module in pvfsga.db.orm and is a dictionary pointing to lists
of dictionaries. The key of the dictionary is the class in the module and the list
of dictionaries contain items to add to that class' table. Each item starts at id
10,000 because the table ids are autoincremented and this is a large enough number to
prevent seeded data from clashing with data added during tests. To expand further,
data added during tests will increment starting at index 1 and we don't want added data
ids to match the seeded data as it will throw an unwanted error. 

We could add the seeded data without an id so the autoincrement starts with this seed
data, but we may want to know for sure what the id is of some data we add in the future.
If we ever expand the seed data, then it may result in needing to change all previously
written tests where the id of the non-seeded data was important.
