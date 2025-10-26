function saveAuthToCookie(value) {
  document.cookie = `scap_auth=${value}; max-age=1800`;
}

function getAuthFromCookie() {
  return document.cookie.replace(
    /(?:(?:^|.*;\s*)scap_auth\s*=\s*([^;]*).*$)|^.*$/,
    "$1"
  );
}

function saveUserInfoToCookie(value) {
  const userInfo = JSON.stringify(value);
  document.cookie = `scap_user=${userInfo}; max-age=1800`;
}

function getUserInfoFromCookie() {
  const userInfo = document.cookie.replace(
    /(?:(?:^|.*;\s*)scap_user\s*=\s*([^;]*).*$)|^.*$/,
    "$1"
  );
  try {
    return userInfo ? JSON.parse(userInfo) : null;
  } catch (e) {
    return null;
  }
}

function deleteCookie() {
  document.cookie = "scap_auth=; expires=Thu, 01 Jan 1970 00:00:01 GMT;";
  document.cookie = "scap_user=; expires=Thu, 01 Jan 1970 00:00:01 GMT;";
}

export { saveAuthToCookie, getAuthFromCookie, deleteCookie, saveUserInfoToCookie, getUserInfoFromCookie };
