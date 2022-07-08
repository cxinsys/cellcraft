function saveAuthToCookie (value) {
  document.cookie = `scap_auth=${value}`
}

function getAuthFromCookie () {
  return document.cookie.replace(
    /(?:(?:^|.*;\s*)scap_auth\s*=\s*([^;]*).*$)|^.*$/,
    '$1'
  )
}

function deleteCookie (value) {
  document.cookie = `${value}=; expires=Thu, 01 Jan 1970 00:00:01 GMT;`
}

export {
  saveAuthToCookie,
  getAuthFromCookie,
  deleteCookie
}
