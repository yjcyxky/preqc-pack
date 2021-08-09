use regex::Regex;

pub fn is_remote_file(remote_path: &str) -> bool {
  let re = Regex::new(r"(oss|s3)://.*$").unwrap();
  re.is_match(remote_path)
}

pub fn parse_remote_path(remote_path: &str) -> (String, String, String) {
  // TODO: how to deal with exception when the filepath is not similar with oss://<bucket-name>/<filepath>
  let re = Regex::new(r"(oss|s3)://([^/]+)(/.*)$").unwrap();
  let cap = re.captures(remote_path).unwrap();
  return (
    String::from(&cap[1]),
    String::from(&cap[2]),
    String::from(&cap[3]),
  );
}