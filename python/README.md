CTS: A Python Implementation
============================

This is a vanilla implementation of the Context Tree Switching algorithm for large alphabets.
Readability is favoured over performance.

Example Usage:

```python
from cts import model

text_model = model.ContextualSequenceModel(context_length=8)

# Open text file for training model. 
with open('alice.txt') as fp:
    # Iterate over text file, byte by byte. 
    while True:
        character = fp.read(1)
        if not character:
            break

        # Train model on this byte. 
        text_model.update(character)

# Now sample 100 more bytes, continuing where the text file ended. 
for t in range(100):
  new_character = text_model.sample()
  # ...
  # Here, do something with sample, e.g. output: 
  # sys.stdout.write(new_character)
  # ...
  # Update observed sequence without training model. 
  text_model.observe(new_character)
```
