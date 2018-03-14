# CWL Tool Instructions

## Installation requirements

* Install `pip` on CentOS 7:
  ```
  sudo yum -y update
  sudo yum -y install python-pip
  ```
* Install `virtualenv` with `pip`:
  ```
  virtualenv -p python2 venv
  ```
* Install `cwltool`, also see [here](https://github.com/common-workflow-language/cwltool):
  ```
  source venv/bin/activate
  pip install cwltool
  ```

## Running CWL Tools
* NOTE:
  * You must fix the file path's in the `.yml` file.
  * `virtualenv` must be active.
    ```
    source venv/bin/activate
    ```

### Run Pilon CWL Tool

```
cwltool pilon.cwl pilon.yml
```
