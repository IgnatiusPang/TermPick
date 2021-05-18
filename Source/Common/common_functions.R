

create_dir_if_not_exists <- function(file_path, mode = "0777") {

  ### Create directory recursively if it doesn't exist
  if (! file.exists(file_path)){

    dir.create(file_path, showWarnings = TRUE, recursive = TRUE, mode = mode)

  }

}


##################################################################################################################

### Function: create_id_to_attribute_hash
### Description: Create a hash function that map keys to attributes.

## Inputs:
## keys: An array of key values
## attributes: An array of attribute values

## Output:
## An environment that act as a hash to convert keys to attributes.

create_id_to_attribute_hash <- function(  keys, attributes) {

  keys <- as.character( as.vector(keys))
  attribute <- as.vector(attributes)

  hash <- new.env(hash = TRUE, parent = parent.frame())

  if ( length(keys) != length(attributes))  {
    warning('Length of keys != Length of attributes list.')
    return(1)
  }

  for ( i in 1:length(keys) )  {
    assign( keys[i], attributes[i], envir = hash)
  }

  return(hash)
}

##################################################################################################################

### Function: convert_key_to_attribute
### Description: Use a predefined hash dictionary to convert any Key to Attribute, return NA if key does not exists

## Inputs:
## key: A key value
## hash: The hash dictionary that maps keys to attributes

## Output:
## A value that correspond to the query key value.

convert_key_to_attribute <- function(key, hash) {

  if ( exists(key, hash) ) {
    return ( get(key, hash))
  }
  else {
    return (NA)
  }
}


##################################################################################################################
